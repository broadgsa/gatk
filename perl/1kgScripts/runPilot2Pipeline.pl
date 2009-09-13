#!/usr/bin/perl -w

use strict;
use Getopt::Long;

sub usage {
    print "Usage: perl runPilot2Pipeline.pl\n\t-i <input dir>\n\t-o <output directory>\n\t[-sting Sting dir]\n\t[-q farm queue; default:gsa]\n\t[-dry]\n";
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

my @samples = ("NA19238","NA19239","NA19240","NA12878","NA12891","NA12892");

foreach my $sample (@samples) {

    my $inputBamSLX = "$inputDir/$sample.pilot2.SLX.bam";
    my $outputHead = "$outputDir/$sample.SLX";
    clean($inputBamSLX, $outputDir, "$outputHead.bam", $queue, $sting, $dry);
#    call($outputBam, $outputHead, $queue, $sting, $dry);

    if ($sample eq "NA12878" || $sample eq "NA19240") {
	my $inputBamSOLID = "$inputDir/$sample.pilot2.SOLID.bam";
	$outputHead = "$outputDir/$sample.SOLID";
	clean($inputBamSOLID, $outputDir, "$outputHead.bam", $queue, $sting, $dry);
#	call($outputBam, $outputHead, $queue, $sting, $dry);

	my $inputBam454 = "$inputDir/$sample.pilot2.454.bam";
	$outputHead = "$outputDir/$sample.454";
#	call($inputBam454, $outputHead, $queue, $sting, $dry);

	my @inputBams = ($inputBamSOLID, $inputBam454);
	$outputHead = "$outputDir/$sample.SOLID_454";
#	multicall($inputBams, $outputHead, $queue, $sting, $dry);

	@inputBams = ($inputBamSLX, $inputBamSOLID, $inputBam454);
	$outputHead = "$outputDir/$sample.allTechs";
#	multicall($inputBams, $outputHead, $queue, $sting, $dry);
    }
}

sub clean {

    my $inputBam = $_[0];
    my $outputDir = $_[1];
    my $outputBam = $_[2];
    my $queue = $_[3];
    my $sting = $_[4];
    my $dry = $_[5];

    my $cmd = "perl $sting/perl/1kgScripts/runCleaningPipeline.pl -i $inputBam -odir $outputDir -obam $outputBam -q $queue -inject -j $outputBam.cleaningpipeline -sting $sting";
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

    my $cmd = "perl $sting/perl/1kgScripts/runCallingPipeline.pl -i $inputBam -o $outputHead -q $queue -snps -wait $inputBam.cleaningpipeline -sting $sting";
    if ($dry) {
	$cmd .= " -dry";
    }
    system($cmd);
}
