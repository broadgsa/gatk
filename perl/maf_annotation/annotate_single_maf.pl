#!/util/bin/perl -w
use strict;

if (scalar(@ARGV) != 1) {
    die("Usage: annotate_single_maf.pl <maf>\n");
}

my ($maf) = @ARGV;

my $cmd = "use matlab && " .
		      "matlab <<STOP\n" .
          "addpath /home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/mike" . "\n" .
          "addpath /home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq" . "\n" .
          "addpath /home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab" . "\n" .
          "annotate_maflite('$maf', '$maf.annotated')" . "\n" .
          "quit" . "\n" .
          "STOP" . "\n";

system($cmd) == 0 or die();
