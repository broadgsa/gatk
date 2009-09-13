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

    $inputBam = "$inputDir/low_coverage_$sample.bam";
    $outputHead = "$outputDir/$sample";
    clean($inputBam, $outputDir, "$outputHead.bam", $queue, $sting, $dry);
#    call($outputBam, $outputHead, $queue, $sting, $dry);
}

sub clean {

    my $inputBam = $_[0];
    my $outputDir = $_[1];
    my $outputBam = $_[2];
    my $queue = $_[3];
    my $sting = $_[4];
    my $dry = $_[5];

    my $cmd = "perl $sting/perl/1kgScripts/runCleaningPipeline.pl -i $inputBam -odir $outputDir -obam $outputBam -q $queue -j $outputBam.cleaningpipeline -sting $sting";
    if ($dry) {
	$cmd .= " -dry";
    }
    system($cmd);
}

sub call {

    my $inputBam = $_[0];
    my $outputHead = $_[1];
    my $queue = $_[2];
    my $sting = $_[3];
    my $dry = $_[4];

    $cmd = "perl $sting/perl/1kgScripts/runCallingPipeline.pl -i $inputBam -o $outputHead -q $queue -wait $inputBam.cleaningpipeline -sting $sting";
    if ($dry) {
	$cmd .= " -dry";
    }
    system($cmd);
}
