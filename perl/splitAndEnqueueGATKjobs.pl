#!/usr/bin/perl -w

# Runs a given GATK walker genome-wide, but splits up the jobs
# and then merges the results together.  One should really use
# scatter-gather here, but I wanted to test against that framework
# (to ensure consistency) so I wrote this temporary script.
# Also, it allowed me to add in farm commands for wait ids
# (which is currently unavailable in scatter-gather).
# If intervals file is left blank, it splits by chromosome.

use strict;
use Getopt::Long;

sub usage {
    print "Usage: perl enqueueGATKjobsByChromosome.pl\n\t-cmd <GATK command without output or -L args>\n\t-o <output file>\n\t[-oarg GATK output argument; default:o]\n\t[-i intervals file]\n\t[-n number of splits to make]\n\t[-q farm queue; default:gsa]\n\t[-wait farm wait id]\n\t[-job farm job name]\n\t[-dry]\n\t[-bam is output format bam?; default:no]\n";
    exit(1);
}

my $cmd = undef;
my $output = undef;
my $outputArg = "o";
my $jobName = undef;
my $wait = undef;
my $queue = "gsa";
my $breaks = 24;
my $intervalsFile = undef;
my $dry;
my $isBam;

GetOptions( "cmd=s" => \$cmd,
	    "o=s" => \$output,
	    "oarg:s" => \$outputArg,
	    "job:s" => \$jobName,
	    "wait:s" => \$wait,
	    "n:i" => \$breaks,
	    "i:s" => \$intervalsFile,
	    "dry!" => \$dry,
	    "bam!" => \$isBam,
	    "q:s" => \$queue);

usage() if ( !$cmd || !$output );

my @intervals;
if (!$intervalsFile) {
    my $num = 1;
    while ($num < 23) {
	push(@intervals, $num);
	$num++;
    }
    push(@intervals, "X");
    push(@intervals, "Y");
} else {
    open(FILE, "< $intervalsFile") or die "can't open $intervalsFile: $!";
    my @lines = <FILE>;
    chomp(@lines);
    my $linecount = scalar(@lines);
    if ($linecount < $breaks) {
	$breaks = $linecount;
    }

    my $linesPerJob = $linecount / $breaks;
    my $index = 0;
    for (my $i = 1; $i < $breaks; $i++) {
	my $interval = "";
	for (my $j = 0; $j < $linesPerJob; $j++) {
	    $interval .= "$lines[$index];";
	    $index++;
	}
	push(@intervals, $interval);
    }
    my $interval = "";
    while ($index < $linecount) {
	$interval .= "$lines[$index];";
	$index++;
    }
    push(@intervals, $interval);
}

my $intervalcount = scalar(@intervals);
for (my $i = 0; $i < $intervalcount; $i++) {
    enqueue($intervals[$i], $cmd, $output, $outputArg, $wait, $queue, $dry, $isBam, ($i+1));
}

mergeResults($output, $queue, $dry, $isBam, $jobName, $intervalcount);

sub enqueue {

    my $interval = $_[0];
    my $cmd = $_[1];
    my $outArg = $_[3];
    my $waitid = $_[4];
    my $queue = $_[5];
    my $dry = $_[6];
    my $bam = $_[7];
    my $index = $_[8];

    my $output = "$_[2].$index";
    if ($bam) {
	$output .= ".bam";
    }

    my $bsub = "bsub -q $queue -o $output.sdout -J $_[2]";
    if ($waitid) {
	$bsub .= " -w \"ended($waitid)\"";
    }

    my $command = "$bsub $cmd -$outputArg $output -L $interval";
    execute($command, $dry);
}

sub mergeResults {

    my $output = $_[0];
    my $queue = $_[1];
    my $dry = $_[2];
    my $bam = $_[3];
    my $jobName = $_[4];
    my $intervalcount = $_[5];

    my $cmd = "bsub -q $queue -o $output.sdout -w \"ended($output)\"";
    if ($jobName) {
	$cmd .= " -J $jobName";
    }

    if ($bam) {
	$cmd .= " samtools merge $output ";
	for (my $i = 1; $i <= $intervalcount; $i++) {
	    $cmd .= "$output.$i.bam ";
	}
    } else {
	$cmd .= " \"cat ";
	for (my $i = 1; $i <= $intervalcount; $i++) {
	    $cmd .= "$output.$i ";
	}
	$cmd .= "> $output\"";
    }

    execute($cmd, $dry);
}

sub execute {

    my $cmd = $_[0];
    my $dry = $_[1];

    if ($dry) {
	print "$cmd\n";
    } else {
	system($cmd);
    }
}
