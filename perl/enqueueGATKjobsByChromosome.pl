#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $walker = undef;
my $pilot = "pilot2";
my $queue = "gsa";
my $tech = "SLX";
my $jar = "/humgen/gsa-scr1/ebanks/Sting/dist/GenomeAnalysisTK.jar";

GetOptions( "T=s" => \$walker,
	    "p:s" => \$pilot,
	    "q:s" => \$queue,
	    "tech:s" => \$tech,
	    "j:s" => \$jar );

exit(1) if ( !$walker );

my @samples;
if ($pilot eq "pilot1") {
    @samples = ("CEU","YRI","CHB-JPT");
} elsif ($pilot eq "pilot2") {
    @samples = ("NA19238","NA19239","NA19240","NA12878","NA12891","NA12892");
}

foreach my $sample (@samples) {

    my $num = 1;
    while ($num < 23) {
	enqueue($sample, $num, $pilot, $queue, $jar, $walker, $tech);
	$num++;
    }

    enqueue($sample, "X", $pilot, $queue, $jar, $walker, $tech);
    enqueue($sample, "Y", $pilot, $queue, $jar, $walker, $tech);
}

sub enqueue {

    my $sample = $_[0];
    my $chr = $_[1];
    my $pilot = $_[2];
    my $queue = $_[3];
    my $jar = $_[4];
    my $walker = $_[5];
    my $tech = $_[6];

    my $inputBam;
    if ($pilot eq "pilot2") {
	$inputBam = "/broad/1KG/DCC/ftp/pilot_data/$sample/alignment/$sample.chrom$chr.$tech.SRP000032.2009_07.bam";
    } else {
	$inputBam = "/broad/hptmp/ebanks/1kg_pilot1/".$sample."_BAMS.list";
    }

    my $outputDir;
    if ($pilot eq "pilot1") {
	$outputDir = "/broad/hptmp/ebanks/1kg_pilot1/cleaned";
    } else {
	$outputDir = "/broad/hptmp/ebanks/1kg_pilot2/cleaned";
    }

    my $cmd;
    my $outputFile;

  SWITCH: {
      $walker eq "IndelIntervals" && do {
	  $outputFile = "$outputDir/intervals/$sample.chr$chr.$tech.indels.intervals";
	  $cmd = "bsub -q $queue -o $outputFile.sdout java -Xmx4096m -jar $jar -S SILENT -T IndelIntervals -R /broad/1KG/reference/human_b36_both.fasta -I $inputBam -o $outputFile -L $chr";
	  last SWITCH;
      };
      $walker eq "MismatchIntervals" && do {
	  $outputFile = "$outputDir/intervals/$sample.chr$chr.$tech.mismatches.intervals";
	  $cmd = "bsub -q $queue -o $outputFile.sdout java -Xmx4096m -jar $jar -S SILENT -T MismatchIntervals -R /broad/1KG/reference/human_b36_both.fasta -I $inputBam -o $outputFile -L $chr";
	  last SWITCH;
      };
      $walker eq "SNPClusters" && do {
	  $outputFile = "$outputDir/intervals/$sample.chr$chr.$tech.clusters.intervals";
	  $cmd = "bsub -q $queue -o $outputFile.sdout java -Xmx4096m -jar $jar -S SILENT -T SNPClusters -R /broad/1KG/reference/human_b36_both.fasta -I $inputBam -o $outputFile -L $chr";
	  last SWITCH;
      };
      $walker eq "IntervalMerger" && do {
	  $outputFile = "$outputDir/intervals/$sample.chr$chr.$tech.merged.intervals";
	  $cmd = "bsub -q $queue -o $outputFile.sdout java -Xmx4096m -jar $jar -S SILENT -T IntervalMerger -R /broad/1KG/reference/human_b36_both.fasta -I $inputBam -o $outputFile -L $chr -intervals $outputDir/intervals/$sample.chr$chr.$tech.indels.intervals -intervals $outputDir/intervals/$sample.chr$chr.$tech.mismatches.intervals";
	  last SWITCH;
      };
      $walker eq "IntervalCleaner" && do {
	  if ($pilot eq "pilot2") {
	      $outputFile = "$outputDir/cleaner/$sample.chr$chr.$tech.bam";
	  } else {
	      $outputFile = "$outputDir/bams/$sample.chr$chr.$tech.bam";
	  }
	  $cmd = "bsub -q $queue -o $outputFile.sdout java -Xmx4096m -jar $jar -S SILENT -T IntervalCleaner -R /broad/1KG/reference/human_b36_both.fasta -I $inputBam -O $outputFile -L $chr -intervals $outputDir/intervals/$sample.chr$chr.$tech.merged.intervals -compress 1";
	  if ($pilot eq "pilot2") {
	      $cmd .= " -cleanedOnly";
	  }
	  last SWITCH;
      };
      $walker eq "CleanedReadInjector" && do {
	  $outputFile = "$outputDir/bams/$sample.chr$chr.$tech.bam";
	  $cmd = "bsub -q $queue -o $outputFile.sdout java -Xmx4096m -jar $jar -S SILENT -T CleanedReadInjector -R /broad/1KG/reference/human_b36_both.fasta -I $inputBam -o $outputFile -L $chr";
	  last SWITCH;
      };

      print "$walker is not a supported class\n";
      exit(1);
    }

#    print "$cmd\n";
    system($cmd);
}
