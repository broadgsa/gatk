#!/usr/bin/perl -w

use strict;
use FindBin;
use lib "$FindBin::Bin/../../";

use StingArgs;
use DistributedMake;
use DataTable;
use Data::Dumper;
use File::Basename;

my %args = &getCommandArguments(
    "LIST" => undef,
    "OUT_DIR" => undef,

    "REFERENCE" => "/broad/1KG/reference/human_b36_both.fasta",
    "DBSNP" => "/humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod",

    "SEQUENCE_INDEX" => "/broad/1KG/DCC/ftp/sequence.index",
    "ALIGNMENT_INDEX" => "/broad/1KG/DCC/ftp/alignment.index",

    "DRY_RUN" => 1,
    "QUEUE" => "gsa",
    "NUM_JOBS" => "",
    "KEEP_GOING" => 1,
);

my $dm = new DistributedMake("dryRun" => $args{'DRY_RUN'}, "queue" => $args{'QUEUE'}, "numJobs" => "", "keepGoing" => 1, "outputFile" => "$args{'OUT_DIR'}/GenerateDataProcessingReport2.log");

my @entries;

open(LIST, $args{'LIST'});
while (my $entry = <LIST>) {
    chomp($entry);

    my %entry;
    ($entry{'POP'}, $entry{'UG'}, $entry{'QCALL'}) = split(/\s+/, $entry);

    push(@entries, \%entry);
}
close(LIST);

my @slides;

foreach my $entry (@entries) {
    my %entry = %{$entry};

    # Combine VCFs
    my $vcfCombineOut = "$args{'OUT_DIR'}/intermediate/$entry{'POP'}.combined.vcf";
    my $vcfCombineOutCmd = "java -Xmx2048m -jar $ENV{'STING_DIR'}/dist/GenomeAnalysisTK.jar -T VCFCombine -R $args{'REFERENCE'} -type UNION -B UG,VCF,$entry{'UG'} -B QCALL,VCF,$entry{'QCALL'} -A -priority 'UG,QCALL' -o $vcfCombineOut -l INFO -L 1";
    $dm->addRule($vcfCombineOut, [$entry{'UG'}, $entry{'QCALL'}], $vcfCombineOutCmd);

    # Evaluate variants
    my $evalRoot = "$vcfCombineOut.eval/eval";
    my $evalCounts = "$evalRoot.Count_Variants.csv";
    my $evalCmd = "java -Xmx2048m -jar $ENV{'STING_DIR'}/dist/GenomeAnalysisTK.jar -T VariantEval -R $args{'REFERENCE'} -D $args{'DBSNP'} -B eval,VCF,$vcfCombineOut -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"UG-filteredInOther\"' -selectName InUG-FilteredInQCALL -select 'set==\"QCALL-filteredInOther\"' -selectName InQCALL-FilteredInUG -select 'set==\"UG\"' -selectName UG -select 'set==\"QCALL\"' -selectName QCALL -reportType R -reportLocation $evalRoot -l INFO -L 1";
    $dm->addRule($evalCounts, $vcfCombineOut, $evalCmd);

    # Get SNPs per sample
    my $snpsPerSample = "$evalRoot.persample.table";
    my $snpsPerSampleCmd = "$ENV{'STING_DIR'}/perl/reports/1KG/GetSNPsPerSample.pl VCF=$vcfCombineOut OUT=$snpsPerSample";
    $dm->addRule($snpsPerSample, $vcfCombineOut, $snpsPerSampleCmd);

    # Create topsheet slide
    my $slideDir = "$args{'OUT_DIR'}/intermediate/$entry{'POP'}/slides";
    my $topSheetRoot = "$slideDir/$entry{'POP'}.000.topsheet";
    my $topSheetSlide = "$topSheetRoot.pdf";
    my $topSheetSlideCmd = "$ENV{'STING_DIR'}/perl/reports/1KG/MakeTopSheetSlide.pl VCF=$vcfCombineOut SEQUENCE_INDEX=$args{'SEQUENCE_INDEX'} ALIGNMENT_INDEX=$args{'ALIGNMENT_INDEX'} EVAL_ROOT=$evalRoot OUT=$topSheetRoot";
    #$dm->addRule($topSheetSlide, [$vcfCombineOut, $args{'SEQUENCE_INDEX'}, $args{'ALIGNMENT_INDEX'}, $evalCounts], $topSheetSlideCmd);

    # Create Venn diagram slide of variants
    my $vennRoot = "$slideDir/$entry{'POP'}.010.venn";
    my $vennPdf = "$vennRoot.pdf";
    my $vennPdfCmd = "$ENV{'STING_DIR'}/perl/reports/1KG/MakeVennSlide.pl EVAL_ROOT=$evalRoot OUT_ROOT=$vennRoot POP=$entry{'POP'}";
    $dm->addRule($vennPdf, $evalCounts, $vennPdfCmd);

    # Create plots for other slides
    my $plotRoot = "$evalRoot.plot";
    my $plotDummy = "$plotRoot.venn_per_ac.norm.pdf";
    my $plotCmd = "Rscript $ENV{'STING_DIR'}/perl/reports/1KG/MakeDataProcessingReportPlots.R $evalRoot $plotRoot $entry{'POP'}";
    $dm->addRule($plotDummy, [$evalCounts, $snpsPerSample], $plotCmd, 'queue' => "");

    # Create other slides
    my $dprSlideDummy = "$slideDir/$entry{'POP'}.050.eval.plot.titv.pdf";
    my $dprSlidesCmd = "$ENV{'STING_DIR'}/perl/reports/1KG/MakeDPRSlides.pl PLOT_ROOT=$plotRoot OUT_DIR=$slideDir POP=$entry{'POP'}";
    $dm->addRule($dprSlideDummy, $plotDummy, $dprSlidesCmd);

    push(@slides, $dprSlideDummy);
}

my $assembledSlides = "$args{'OUT_DIR'}/dpr.pdf";
my $assembledSlidesCmd = "$ENV{'STING_DIR'}/perl/reports/1KG/AssembleSlides.pl DIR=$args{'OUT_DIR'} OUT=$assembledSlides";
$dm->addRule($assembledSlides, \@slides, $assembledSlidesCmd);

$dm->execute();
