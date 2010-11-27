package DistributedMake;

use strict;
use File::Temp qw/ tempfile tempdir /;
use File::Basename;

sub parseHostsString {
    my ($hoststring) = @_;

    if ($hoststring !~ /\s+\+\s+/) {
        return undef;
    }

    my @hostobjs = split(/\s+\+\s+/, $hoststring);

    my @hosts;
    foreach my $hostobj (@hostobjs) {
        my ($multiplier, $server) = $hostobj =~ /(\d+)\*(\w+)/;
        for (my $i = 0; $i < $multiplier; $i++) {
            push(@hosts, $server);
        }
    }

    return \@hosts;
}

sub new {
    my ($class, %args) = @_;

    my %self = (
        'dryRun'             => 1,
        'numJobs'            => undef,
        'keepGoing'          => 0,
        'alwaysMake'         => 0,
        'debugging'          => 0,
        'ignoreErrors'       => 0,
        'printDirectory'     => 0,
        'unlink'             => 1,
        'hosts'              => "",

        'queue'              => undef,
        'memLimit'           => 2,
        'outputFile'         => 'distributedmake.log',
        'mailTo'             => 'crd-lsf@broad.mit.edu',
        'wait'               => 1,
        'rerunnable'         => 0,
        'migrationThreshold' => undef,
        'extra'              => '',

        'target'             => 'all',

        %args,

        'targets'            => [],
        'hostindex'          => 0,
    );

    $self{'makefile'}  = new File::Temp(TEMPLATE => "/tmp/DistributedMake_XXXXXX", SUFFIX => ".makefile", UNLINK => $self{'unlink'}),
    $self{'hostarray'} = &parseHostsString($self{'hosts'});
    $self{'projectName'}   = basename($self{'makefile'});

	bless \%self, $class;

	return \%self;
}

sub addRule {
    my ($self, $targetsref, $dependenciesref, $cmdsref, %batchjoboverrides) = @_;
    my @targets      = (ref($targetsref)      eq 'ARRAY') ? @$targetsref      : ( $targetsref );
    my @dependencies = (ref($dependenciesref) eq 'ARRAY') ? @$dependenciesref : ( $dependenciesref );
    my @cmds         = (ref($cmdsref)         eq 'ARRAY') ? @$cmdsref         : ( $cmdsref );

    my @prepcmds;

    my $cmdprefix = "";
    if (defined($self->{'hostarray'})) {
        $cmdprefix = "ssh ${$self->{'hostarray'}}[$self->{'hostindex'}] ";

        $self->{'hostindex'}++;
        if ($self->{'hostindex'} == scalar(@{$self->{'hostarray'}}) - 1) {
            $self->{'hostindex'} = 0;
        }
    } elsif ((defined($self->{'queue'}) && $self->{'queue'} ne '' && (exists($batchjoboverrides{'queue'}) ? defined($batchjoboverrides{'queue'}) && $batchjoboverrides{'queue'} ne '' : 1)) || (exists($batchjoboverrides{'queue'}) && defined($batchjoboverrides{'queue'}) && $batchjoboverrides{'queue'} ne '')) {
        my %bja = (
            'queue' => $self->{'queue'},
            'memLimit' => $self->{'memLimit'},
            'projectName' => $self->{'projectName'},
            'outputFile' => $self->{'outputFile'},
            'mailTo' => $self->{'mailTo'},
            'wait' => $self->{'wait'},
            'rerunnable' => $self->{'rerunnable'},
            'migrationThreshold' => $self->{'migrationThreshold'},
            'extra' => $self->{'extra'},
            %batchjoboverrides,
        );

        my $rerunnable = $bja{'rerunnable'} ? "-r" : "";
        my $migrationThreshold = $bja{'rerunnable'} && defined($bja{'migrationThreshold'}) ? "-mig $bja{'migrationThreshold'}" : "";
        my $wait = $bja{'wait'} ? "-K" : "";

        my $logdir = dirname($bja{'outputFile'});
        if (!-e $logdir) {
            my $mklogdircmd = "\@test \"! -d $logdir\" && mkdir -p $logdir";
            push(@prepcmds, $mklogdircmd);
        }

        my $memRequest = $bja{'memLimit'} * 1.5;
        my $integerMemRequest = int($memRequest);
        my $memCutoff = $bja{'memLimit'} * 1024 * 1024 * 1.25;

        # A quick check to make sure that java commands being dispatched to the farm are instructed to run under a default memory limit
        for (my $i = 0; $i <= $#cmds; $i++) {
            if ($cmds[$i] =~ /^java / && $cmds[$i] =~ / -jar / && $cmds[$i] !~ / -Xmx/) {
                $cmds[$i] =~ s/^java /java -Xmx$bja{'memLimit'}g /;
            }
        }

        $cmdprefix = "bsub -q $bja{'queue'} -M $memCutoff -P $bja{'projectName'} -o $bja{'outputFile'} -u $bja{'mailTo'} -R \"rusage[mem=$integerMemRequest]\" $wait $rerunnable $migrationThreshold $bja{'extra'}    ";
    }

    my $rootdir = dirname($targets[0]);
    if (!-e $rootdir) {
        my $mkdircmd = "\@test \"! -d $rootdir\" && mkdir -p $rootdir";
        push(@prepcmds, $mkdircmd);
    }

    # We have to touch the final file just in case the time between different nodes on the farm are not synchronized.
    print { $self->{'makefile'} } "$targets[0]: " . join(" ", @dependencies) . "\n\t" . join("\n\t", @prepcmds) . "\n\t$cmdprefix" . join("\n\t$cmdprefix", @cmds) . "\n\ttouch -c $targets[0]\n\n\n";

    push(@{$self->{'targets'}}, $targets[0]);
}

sub execute {
    my ($self, %overrides) = @_;

    print { $self->{'makefile'} } "all: " . join(" ", @{$self->{'targets'}}) . "\n\n";
    print { $self->{'makefile'} } ".DELETE_ON_ERROR:\n";

	my %makeargs = (
        'dryRun'         => $self->{'dryRun'},
        'numJobs'        => $self->{'numJobs'},
        'keepGoing'      => $self->{'keepGoing'},
        'alwaysMake'     => $self->{'alwaysMake'},
        'debugging'      => $self->{'debugging'},
        'ignoreErrors'   => $self->{'ignoreErrors'},
        'printDirectory' => $self->{'printDirectory'},
        'target'         => $self->{'target'},
		%overrides,
	);

    my $numjobs = $makeargs{'numJobs'};
    if (!defined($numjobs)) {
        if (defined($self->{'hostarray'}) && scalar($self->{'hostarray'}) > 0) {
            $numjobs = scalar(@{$self->{'hostarray'}});
        } else {
            $numjobs = 1;
        }
    }

    my $makecmd = "make" .
                    ($makeargs{'dryRun'}         ? " -n" : "") .
                    ($makeargs{'keepGoing'}      ? " -k" : "") .
                    ($makeargs{'alwaysMake'}     ? " -B" : "") .
                    ($makeargs{'ignoreErrors'}   ? " -i" : "") .
                    ($makeargs{'printDirectory'} ? " -w" : "") .
                    ($makeargs{'debugging'} =~ /[abvijm]+/ ? " --debug=$makeargs{'debugging'}" : "") .
                    ($makeargs{'debugging'} =~ /\d+/ && $makeargs{'debugging'} == 1 ? " -d" : "") .
                    " -j $numjobs" .
                    " -f " . $self->{'makefile'}->filename .
                    " $makeargs{'target'}";

    print "$makecmd\n";
    system($makecmd);
    print "$makecmd\n";
}

1;
