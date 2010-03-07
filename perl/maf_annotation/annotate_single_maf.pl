#!/util/bin/perl -w
use strict;

if (scalar(@ARGV) != 1) {
    die("Usage: annotate_single_maf.pl <maf>\n");
}

my ($maf) = @ARGV;

# TODO: Have a common checkout of https://svn/CancerGenomeAnalysis
# or a compiled version of this matlab program like Firehose uses.
my $cmd = "matlab <<STOP\n" .
          "addpath /humgen/gsa-hphome1/kshakir/src/cga_matlab_r7169/mike" . "\n" .
          "addpath /humgen/gsa-hphome1/kshakir/src/cga_matlab_r7169/seq" . "\n" .
          "addpath /humgen/gsa-hphome1/kshakir/src/cga_matlab_r7169" . "\n" .
          "annotate_maflite('$maf', '$maf.annotated')" . "\n" .
          "quit" . "\n" .
          "STOP" . "\n";

system($cmd) == 0 or die();
