#!/usr/bin/perl -w

use strict;
use Getopt::Long;

sub usage {
    print "Usage: perl runCleaningPipeline.pl\n\t-i <input bam>\n\t-obam <output bam name>\n\t[-sting Sting dir]\n\t[-badsnps badsnps file name]\n\t[-inject]\n\t[-q farm queue; default:gsa]\n\t[-wait farm wait id]\n\t[-job final farm job name]\n\t[-dry]\n";
    exit(1);
}

my $inputBam = undef;
my $outputBam = undef;
my $jobName = undef;
my $wait = undef;
my $badsnps = undef;
my $dry;
my $inject;
my $queue = "gsa";
my $sting = "/humgen/gsa-scr1/ebanks/Sting";

GetOptions( "i=s" => \$inputBam,
	    "obam=s" => \$outputBam,
	    "q:s" => \$queue,
	    "dry!" => \$dry,
	    "inject!" => \$inject,
	    "job:s" => \$jobName,
	    "wait:s" => \$wait,
	    "badsnps:s" => \$badsnps,
	    "sting:s" => \$sting );

usage() if ( !$inputBam || !$outputBam );

my $indelIntervals = "$outputBam.indels.intervals";
my $command = "perl $sting/perl/splitAndEnqueueGATKjobs.pl -cmd \"java -Xmx4096m -jar $sting/dist/GenomeAnalysisTK.jar -S SILENT -T IndelIntervals -R /broad/1KG/reference/human_b36_both.fasta -I $inputBam\" -o $indelIntervals -oarg o -j $outputBam.intervals -q $queue";
if ($dry) {
    $command .= " -dry";
}
if ($wait) {
    $command .= " -wait $wait";
}
system($command);

my $mismatchIntervals = "$outputBam.mismatches.intervals";
$command = "perl $sting/perl/splitAndEnqueueGATKjobs.pl -cmd \"java -Xmx4096m -jar $sting/dist/GenomeAnalysisTK.jar -S SILENT -T MismatchIntervals -R /broad/1KG/reference/human_b36_both.fasta -I $inputBam\" -o $mismatchIntervals -oarg o -j $outputBam.intervals -q $queue";
if ($dry) {
    $command .= " -dry";
}
if ($wait) {
    $command .= " -wait $wait";
}
system($command);

my $mergedIntervals = "$outputBam.merged.intervals";
$command = "perl $sting/perl/splitAndEnqueueGATKjobs.pl -cmd \"java -Xmx4096m -jar $sting/dist/GenomeAnalysisTK.jar -S SILENT -T IntervalMerger -R /broad/1KG/reference/human_b36_both.fasta -I $inputBam -intervals $indelIntervals -intervals $mismatchIntervals\" -o $mergedIntervals -oarg o -wait $outputBam.intervals -j $outputBam.merged -q $queue";
if ($dry) {
    $command .= " -dry";
}
system($command);

my $cleanedBam;
if ($inject) {
    $cleanedBam = "$outputBam.cleaned.bam";
} else {
    $cleanedBam = "$outputBam";
}
$command = "perl $sting/perl/splitAndEnqueueGATKjobs.pl -cmd \"java -Xmx4096m -jar $sting/dist/GenomeAnalysisTK.jar -S SILENT -T IntervalCleaner -R /broad/1KG/reference/human_b36_both.fasta -I $inputBam -compress 1";
if ($inject) {
    $command .= " -cleanedOnly\" -j $outputBam.cleaned";
} else {
    $command .= "\" -j $jobName";
}
$command .= " -o $cleanedBam -oarg O -wait $outputBam.merged -q $queue -bam -i $mergedIntervals -n 50";
if ($dry) {
    $command .= " -dry";
}
system($command);

my $snpsFile = $badsnps;
if (!$snpsFile) {
    $snpsFile = "$outputBam.badsnps";
}
$command = "perl $sting/perl/splitAndEnqueueGATKjobs.pl -cmd \"java -Xmx4096m -jar $sting/dist/GenomeAnalysisTK.jar -S SILENT -T IntervalCleaner -R /broad/1KG/reference/human_b36_both.fasta -I $inputBam\"";
if ($inject) {
    $command .= " -j $outputBam.cleaned";
} else {
    $command .= " -j $jobName";
}
$command .= "\" -o $snpsFile -oarg snps -wait $outputBam.merged -q $queue -i $mergedIntervals -n 50";
if ($dry) {
    $command .= " -dry";
}
system($command);

if ($inject) {
    my $bam = "$outputBam";
    $command = "bsub -q $queue -o $bam.sdout";
    if ($jobName) {
	$command .= " -J $jobName";
    }
    $command .= " -w \"ended($outputBam.cleaned)\" java -Xmx4096m -jar $sting/dist/GenomeAnalysisTK.jar -S SILENT -T CleanedReadInjector -R /broad/1KG/reference/human_b36_both.fasta -I $inputBam --output_bam $bam --cleaned_reads $cleanedBam -compress 1";
    if ($dry) {
	print "$command\n";
    } else {
	system($command);
    }
}
