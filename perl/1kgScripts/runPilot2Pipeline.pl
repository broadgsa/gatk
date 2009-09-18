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
	    "o=s" => \$outputDir,
	    "q:s" => \$queue,
	    "dry!" => \$dry,
	    "sting:s" => \$sting );

usage() if ( !$inputDir || !$outputDir );

my @samples = ("NA19238","NA19239","NA19240","NA12878","NA12891","NA12892");

# Official genome-wide Depth of Coverage tables for pilot 2, freeze 5:
#        NA12878 NA12891 NA12892 NA19238 NA19239 NA19240
# 454:     36                                      18
# SLX:     82      91      70      56      68      86
# SOLID:   37                                      64
# 454+SLD: 64                                      77
# ALL:     138                                     150

my %DoC_454 = ( "NA12878" => 36, "NA19240" => 18 );
my %DoC_slx = ( "NA12878" => 82, "NA12891" => 91, "NA12892" => 70,"NA19238" => 56, "NA19239" => 68, "NA19240" => 86 );
my %DoC_solid = ( "NA12878" => 37, "NA19240" => 64 );
my %DoC_454solid = ( "NA12878" => 64, "NA19240" => 77 );
my %DoC_all = ( "NA12878" => 138, "NA19240" => 150 );
my %MQ_hash = ( "SLX" => 100, "SOLID" => 5, "454" => 5, "454SOLID" => 10, "ALL" => 110 );

foreach my $sample (@samples) {

    my $inputBamSLX = "$inputDir/$sample.pilot2.SLX.bam";
    my $outputHeadSLX = "$outputDir/$sample.SLX";
    my $outputBamSLX = "$outputHeadSLX.bam";
    my $badsnpsSLX = "$outputBamSLX.badsnps";
    clean($inputBamSLX, $outputBamSLX, $queue, $sting, $dry, $badsnpsSLX);
    call("-I $outputBamSLX", $outputHeadSLX, $queue, $sting, $dry, "$inputBamSLX.cleaningpipeline", $sample, $badsnpsSLX, $DoC_slx{$sample}, $MQ_hash{"SLX"});

    if ($sample eq "NA12878" || $sample eq "NA19240") {
	my $inputBamSOLID = "$inputDir/$sample.pilot2.SOLID.bam";
	my $outputHeadSOLID = "$outputDir/$sample.SOLID";
	my $outputBamSOLID = "$outputHeadSOLID.bam";
	my $badsnpsSOLID = "$outputBamSOLID.badsnps";
	clean($inputBamSOLID, $outputBamSOLID, $queue, $sting, $dry, $badsnpsSOLID);
	call("-I $outputBamSOLID", $outputHeadSOLID, $queue, $sting, $dry, "$inputBamSOLID.cleaningpipeline", $sample, $badsnpsSOLID, $DoC_solid{$sample}, $MQ_hash{"SOLID"});

	my $inputBam454 = "$inputDir/$sample.pilot2.454.bam";
	my $outputHead454 = "$outputDir/$sample.454";
	call("-I $inputBam454", $outputHead454, $queue, $sting, $dry, undef, $sample, $badsnpsSLX, $DoC_454{$sample}, $MQ_hash{"454"});

	my $outputHead = "$outputDir/$sample.SOLID_454";
	call("-I $outputBamSOLID -I $inputBam454", $outputHead, $queue, $sting, $dry, "$inputBamSOLID.cleaningpipeline", $sample, $badsnpsSOLID, $DoC_454solid{$sample}, $MQ_hash{"454SOLID"});

	$outputHead = "$outputDir/$sample.allTechs";
	call("-I $outputBamSLX -I $outputBamSOLID -I $inputBam454", $outputHead, $queue, $sting, $dry, "*.cleaningpipeline", $sample, $badsnpsSLX, $DoC_all{$sample}, $MQ_hash{"ALL"});
    }
}

sub clean {

    my $inputBam = $_[0];
    my $outputBam = $_[1];
    my $queue = $_[2];
    my $sting = $_[3];
    my $dry = $_[4];
    my $badsnps = $_[5];

    my $cmd = "perl $sting/perl/1kgScripts/runCleaningPipeline.pl -i $inputBam -obam $outputBam -q $queue -inject -j $outputBam.cleaningpipeline -sting $sting -badsnps $badsnps";
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
    my $doc = $_[8];
    my $mq = $_[9];

    my $cmd = "perl $sting/perl/1kgScripts/runCallingPipeline.pl -i \"$inputBams\" -o $outputHead -q $queue -snps -sting $sting -sample $sample -badsnps $badsnps -doc $doc -mq $mq";
    if ($dry) {
	$cmd .= " -dry";
    }
    if ($wait) {
	$cmd .= " -wait $wait";
    }
    system($cmd);
}
