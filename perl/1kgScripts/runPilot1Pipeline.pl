#!/usr/bin/perl -w

use strict;
use Getopt::Long;

sub usage {
    print "Usage: perl runPilot1Pipeline.pl\n\t-i <input dir>\n\t-o <output directory>\n\t[-sting Sting dir]\n\t[-q farm queue; default:gsa]\n\t[-dry]\n";
    exit(1);
}

my $inputDir = undef;
my $outputDir = undef;
my $dry;
my $queue = "gsa";
my $sting = "/humgen/gsa-scr1/ebanks/Sting";

GetOptions( "i=s" => \$inputDir,
	    "odir=s" => \$outputDir,
	    "q:s" => \$queue,
	    "dry!" => \$dry,
	    "sting:s" => \$sting );

usage() if ( !$inputDir || !$outputDir );

my @samples = ("ceu","yri","chb_jpt");

foreach my $sample (@samples) {

    my $inputBam = "$inputDir/low_coverage_$sample.bam";
    my $outputHead = "$outputDir/$sample";
    my $outputBam = "$outputHead.bam";
    my $badsnps = "$outputBam.badsnps";
    clean($inputBam, $outputBam, $queue, $sting, $dry, $badsnps);
    call("-I $outputBam", $outputHead, $queue, $sting, $dry, "$inputBam.cleaner.*", $sample, $badsnps);
}

sub clean {

    my $inputBam = $_[0];
    my $outputBam = $_[1];
    my $queue = $_[2];
    my $sting = $_[3];
    my $dry = $_[4];
    my $badsnps = $_[5];

    my $cmd = "perl $sting/perl/1kgScripts/runCleaningPipeline.pl -i $inputBam -obam $outputBam -q $queue -j $outputBam.cleaner.pipeline -sting $sting -badsnps $badsnps";
    if ($dry) {
	$cmd .= " -dry";
    }
    system($cmd);
}

sub call {

    my $inputBams = $_[0];
    my $outputHead = $_[1];
    my $queue = $_[2];
    my $sting = $_[3];
    my $dry = $_[4];
    my $wait = $_[5];
    my $sample = $_[6];
    my $badsnps = $_[7];

    my $cmd = "perl $sting/perl/1kgScripts/runCallingPipeline.pl -i $inputBams -o $outputHead -q $queue -sting $sting -frac 0 -sample $sample -badsnps $badsnps";
    if ($dry) {
	$cmd .= " -dry";
    }
    if ($wait) {
	$cmd .= " -wait $wait";
    }
    system($cmd);
}
