#!/usr/bin/perl -w

use strict;
use lib "$ENV{'STING_DIR'}/perl";
use StingArgs;
use File::Basename;
use DistributedMake;

my %args = &getCommandArguments("PLOT_ROOT" => undef, "OUT_DIR" => undef, "POP" => undef);

my $dm = new DistributedMake("dryRun" => 0, "numJobs" => "1");

my $timestamp = gmtime();

my @items = (
    { 'plotprefix' => "$args{'PLOT_ROOT'}.venn_per_ac", 'title' => 'Callset concordance per allele count', 'population' => $args{'POP'}, 'caption' => "Concordance of the UnifiedGenotyper and QCALL callsets, per allele count.  Bottom and top sections of each bar represent the calls unique to each caller's callset, while the middle portion represents the intersection.", 'timestamp' => $timestamp },
    { 'plotprefix' => "$args{'PLOT_ROOT'}.venn_per_ac.norm", 'title' => 'Normalized callset concordance per allele count', 'population' => $args{'POP'}, 'caption' => "Concordance of the UnifiedGenotyper and QCALL callsets, per allele count, normalized by the number of variants in the callset union for that allele count bin.  Bottom and top sections of each bar represent the calls unique to each caller's callset, while the middle portion represents the intersection.", 'timestamp' => $timestamp },
    { 'plotprefix' => "$args{'PLOT_ROOT'}.acs", 'title' => 'Allele count spectrum', 'population' => $args{'POP'}, 'caption' => 'Allele count spectrum for all known, novel, and all variants.  Individual contributions from separate callers (thin, dotted or dashed lines) as well as the callset intersection (bold, solid lines) are shown.', 'timestamp' => $timestamp },
    { 'plotprefix' => "$args{'PLOT_ROOT'}.titv", 'title' => 'Ti/Tv per allele count', 'population' => $args{'POP'}, 'caption' => 'Ratio of transition mutations to transversions, per allele count.  Whole-genome expectation is approximately 2.1.  Fluctuation at higher allele counts is typically due to fewer available variants for the calculation at those frequencies.', 'timestamp' => $timestamp },
    { 'plotprefix' => "$args{'PLOT_ROOT'}.counts_per_sample", 'title' => 'SNP calls per sample', 'population' => $args{'POP'}, 'caption' => 'Number of SNPs per sample.  As samples in plot are all from the same or similar populations, with presumably similar heterozygosity, all samples should have roughly consistent performance.  Outliers may indicate problematic samples - perhaps in coverage, identity QC (sample mixup), or quality of sequencing.', 'timestamp' => $timestamp },
);

my $idnum = 1;

foreach my $item (@items) {
    $idnum++;
    my $id = "$args{'POP'}." . sprintf("%03d", $idnum*10);

    my %item = %$item;

    my $basename = "$id." . basename($item{'plotprefix'});
    my $table = "$args{'OUT_DIR'}/$basename.kvp";

    open(TABLE, ">$table");
    foreach my $key (sort { $a cmp $b } keys(%item)) {
        print TABLE "$key\t$item{$key}\n";
    }
    close(TABLE);

    my $render = "$args{'OUT_DIR'}/$basename.tex";
    my $renderCmd = "$ENV{'STING_DIR'}/perl/reports/1KG/RenderTemplate.pl KVP=$table TEMPLATE=$ENV{'STING_DIR'}/perl/reports/1KG/plot_and_text_template.tex TEMPLATE_OUT=$render";
    $dm->addRule($render, "$item{'plotprefix'}.pdf", $renderCmd);

    my $pdf = "$args{'OUT_DIR'}/$basename.pdf";
    my $pdfdir = dirname($pdf);
    my $pdfCmd = "pdflatex -output-directory $pdfdir $render";
    $dm->addRule($pdf, $render, $pdfCmd);
}

$dm->execute();
