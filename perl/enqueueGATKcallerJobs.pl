#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $pilot = "pilot2";
my $queue = "gsa";
my $tech = "SLX";
my $jar = "/humgen/gsa-scr1/ebanks/Sting/dist/GenomeAnalysisTK.jar";

GetOptions( "p:s" => \$pilot,
	    "q:s" => \$queue,
	    "tech:s" => \$tech,
	    "j:s" => \$jar );

my @samples;
if ($pilot eq "pilot1") {
    @samples = ("CEU","YRI","CHB-JPT");
} elsif ($pilot eq "pilot2") {
    @samples = ("NA19238","NA19239","NA19240","NA12878","NA12891","NA12892");
}

foreach my $sample (@samples) {
    enqueue($sample, $pilot, $queue, $jar, $tech);
}

sub enqueue {

    my $sample = $_[0];
    my $pilot = $_[1];
    my $queue = $_[2];
    my $jar = $_[3];
    my $tech = $_[4];

    my $inputBamStr = "";
    my $outputDir;
    if ($pilot eq "pilot2") {
	$inputBamStr = "-I /humgen/gsa-hphome1/projects/1kg_pilot2/useTheseBamsForAnalyses/$sample.$tech.bam";
	$outputDir = "/broad/hptmp/ebanks/1kg_pilot2/cleaned/calls";
    } else {
	my $num = 1;
	while ($num < 23) {
	    $inputBamStr .= "-I /broad/hptmp/ebanks/1kg_pilot1/cleaned/bams/$sample.chr$num.$tech.bam ";
	    $num++;
	}
	$inputBamStr .= "-I /broad/hptmp/ebanks/1kg_pilot1/cleaned/bams/$sample.chrX.$tech.bam -I /broad/hptmp/ebanks/1kg_pilot1/cleaned/bams/$sample.chrY.$tech.bam ";
	$outputDir = "/broad/hptmp/ebanks/1kg_pilot1/cleaned/calls";
    }

    my $outputFile = "$outputDir/indels/$sample.$tech.low.calls";
    my $cmd = "bsub -q $queue -o $outputFile.sdout java -Xmx4096m -jar $jar -S SILENT -T IndelGenotyper -R /broad/1KG/reference/human_b36_both.fasta $inputBamStr -o $outputFile -minConsensusFraction 0.5 -minFraction 0.";
    if ($pilot eq "pilot1") { $cmd .= "0"; }
    $cmd .= "1 -minCnt 2 -1kg";
    system($cmd);

    $outputFile = "$outputDir/indels/$sample.$tech.high.calls";
    $cmd = "bsub -q $queue -o $outputFile.sdout java -Xmx4096m -jar $jar -S SILENT -T IndelGenotyper -R /broad/1KG/reference/human_b36_both.fasta $inputBamStr -o $outputFile -minConsensusFraction 0.5 -minFraction 0.";
    if ($pilot eq "pilot1") { $cmd .= "0"; }
    $cmd .= "3 -minCnt 2 -1kg";
    system($cmd);

    if ($pilot eq "pilot2") {
	$outputFile = "$outputDir/unfiltered_snps/$sample.$tech.geli.calls";
	$cmd = "bsub -q $queue -o $outputFile.sdout java -Xmx4096m -jar $jar -S SILENT -T SingleSampleGenotyper -R /broad/1KG/reference/human_b36_both.fasta $inputBamStr -varout $outputFile -lod 5";
	system($cmd);
    }
}
