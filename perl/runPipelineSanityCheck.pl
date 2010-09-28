#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Data::Dumper;
use File::Path;

my @log;

sub getTimestamp {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    my $timestamp = sprintf("%04d.%02d.%02d", 1900 + $year, $mon, $mday);

    return $timestamp;
}

sub log {
    my ($message) = @_;

    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    my $timestamp = sprintf("[%04d.%02d.%02d %02d:%02d:%02d]", 1900 + $year, $mon, $mday, $hour, $min, $sec);

    my @message = split(/\n/, $message);
    for (my $i = 0; $i <= $#message; $i++) {
        $message[$i] = "$timestamp $message[$i]";
    }

    print STDOUT join("\n", @message) . "\n";

    push(@log, @message);
}

sub emailLog {
    my ($subject, $logfile, $email) = @_;

    my $cmd = "mail -s 'Cron <kiran\@gsa3> $subject' $email < $logfile";
    system($cmd);
}

sub killPipeline {
    chomp(my @jobs = qx(bjobs -P QueuePipelineTest));

    my @ids;
    foreach my $job (@jobs) {
        if ($job =~ /^(\d+)\s+/) {
            push(@ids, $1);
        }
    }

    system("bkill " . join(" ", @ids));
}

sub execute {
    my ($command, $output, $timeout, $run, $logfile, $email) = @_;

    my $status = 0;
    my $result = "(output not generated in dry-run mode)";

    &log("Execute: '$command'");

    if ($run == 1) {
        eval {
            local $SIG{'ALRM'} = sub { die("alarm\n") };
            alarm($timeout);

            $result = qx($command 2>&1);
            $status = $?;

            alarm(0);
        };

        if ($@) {
            &log("Timeout: '$command' did not return within '$timeout' seconds");

            if ($run == 1) {
                open(LOG, ">$logfile");
                print LOG join("\n", @log);
                close(LOG);

                &emailLog("Analysis pipeline nightly test: failed", $logfile, $email);
                &killPipeline();
            }

            exit($status);
        }
    }

    if ($status == 0) {
        if ($output == 0) {
            &log("Success: '$command' exited with code '$status'");
        } else {
            &log("Success: '$command' exited with code '$status'  Output:");
            &log($result);
        }

        return 0;
    }

    &log("Failure: '$command' exited with code '$status'  Output:");

    if (!defined($result)) {
        $result = "(unable to capture)\n";
    }

    &log($result);

    if ($run == 1) {
        open(LOG, ">$logfile");
        print LOG join("\n", @log);
        close(LOG);

        &emailLog("Analysis pipeline nightly test: failed", $logfile, $email);
        &killPipeline();
    }

    exit($status);
}

# Default arguments
my %args = (
    'pipelineRoot' => "/humgen/gsa-hpprojects/analysis/pipeline",
    'outputRoot' => "/humgen/gsa-hpprojects/analysis/pipeline/projects/QueuePipelineNightlyTest",
    'testRepositoryDir' => "/humgen/gsa-hpprojects/analysis/pipeline/repositories/StingTest",
    'testInputYAML' => "/humgen/gsa-hpprojects/GATK/data/Validation_Data/QueuePipelineTestData/QueuePipelineTestData.yaml",
    'email' => "$ENV{'USER'}\@broadinstitute.org",
    'run' => 0,
);

# Get command-line options
GetOptions(
    'outputRoot=s' => \$args{'outputRoot'},
    'pipelineRoot=s' => \$args{'pipelineRoot'},
    'testRepositoryDir=s' => \$args{'testRepositoryDir'},
    'testInputYAML=s' => \$args{'testInputYAML'},
    'email=s' => \$args{'email'},
    'run' => \$args{'run'},
    'help|h' => \$args{'help'}
);

# Should we emit some help?
if ($args{'help'}) {
    print Dumper(\%args) . "\n";
    exit(0);
}

# Set some directories
$args{'testDir'} = "$args{'outputRoot'}/" . &getTimestamp();
$args{'testOutputDir'} = "$args{'testDir'}/output";
$args{'logFile'} = "$args{'testDir'}/output.log";

# Set up environment
$ENV{'PATH'} = "$args{'testRepositoryDir'}/shell:$args{'testRepositoryDir'}/python:$ENV{'PATH'}";

# Execute jobs
&execute("echo \$PATH",                     1, 10, $args{'run'}, $args{'logFile'}, $args{'email'});
&execute("mkdir -p $args{'testDir'}",       0, 10, $args{'run'}, $args{'logFile'}, $args{'email'});
&execute("mkdir -p $args{'testOutputDir'}", 0, 10, $args{'run'}, $args{'logFile'}, $args{'email'});
&execute("rm -rf $args{'logFile'}",         0, 10, $args{'run'}, $args{'logFile'}, $args{'email'});

&execute("cd $args{'testRepositoryDir'} && svn cleanup",                        0, 120,  $args{'run'}, $args{'logFile'}, $args{'email'});
&execute("cd $args{'testRepositoryDir'} && svn up",                             0, 600,  $args{'run'}, $args{'logFile'}, $args{'email'});
&execute("cd $args{'testRepositoryDir'} && ant clean",                          0, 600,  $args{'run'}, $args{'logFile'}, $args{'email'});
&execute("cd $args{'testRepositoryDir'} && ant dist playground oneoffs queue",  0, 1200, $args{'run'}, $args{'logFile'}, $args{'email'});
&execute("cd $args{'testRepositoryDir'} && svn info",                           1, 10,   $args{'run'}, $args{'logFile'}, $args{'email'});

&execute("cd $args{'testOutputDir'} && java -jar $args{'testRepositoryDir'}/dist/Queue.jar -jobProject QueuePipelineTest -S $args{'testRepositoryDir'}/scala/qscript/fullCallingPipeline.q -Y $args{'testInputYAML'} -refseqTable /humgen/gsa-hpprojects/GATK/data/Annotations/refseq/refGene-big-table-hg18.txt --gatkjar $args{'testRepositoryDir'}/dist/GenomeAnalysisTK.jar -contigIntervals /humgen/gsa-hpprojects/analysis/pipeline/resources/whole_genome.intervals -titv 3.0 -skipCleaning -bsub -bsubWait", 0, 120, $args{'run'}, $args{'logFile'}, $args{'email'});
&execute("cd $args{'testOutputDir'} && java -jar $args{'testRepositoryDir'}/dist/Queue.jar -jobProject QueuePipelineTest -S $args{'testRepositoryDir'}/scala/qscript/fullCallingPipeline.q -Y $args{'testInputYAML'} -refseqTable /humgen/gsa-hpprojects/GATK/data/Annotations/refseq/refGene-big-table-hg18.txt --gatkjar $args{'testRepositoryDir'}/dist/GenomeAnalysisTK.jar -contigIntervals /humgen/gsa-hpprojects/analysis/pipeline/resources/whole_genome.intervals -titv 3.0 -skipCleaning -bsub -bsubWait -run", 0, 1200, $args{'run'}, $args{'logFile'}, $args{'email'});

&log("All tests completed successfully");

if ($args{'run'} == 1) {
    open(LOG, ">$args{'logFile'}");
    print LOG join("\n", @log);
    close(LOG);

    &emailLog("Analysis pipeline nightly test: succeeded", $args{'logFile'}, $args{'email'});
}
