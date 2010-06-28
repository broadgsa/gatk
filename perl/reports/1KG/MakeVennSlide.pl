#!/usr/bin/perl -w

use strict;
use lib "$ENV{'STING_DIR'}/perl";
use StingArgs;
use DistributedMake;
use Data::Dumper;
use File::Basename;

sub getCSV {
    my ($file) = @_;

    my @header;
    my @csv;

    open(CSV, $file);
    while (my $csvline = <CSV>) {
        if ($csvline !~ /#/) {
            chomp(my @columns = split(/,/, $csvline));

            if (scalar(@columns) > 1) {
                if (scalar(@header) < 2) {
                    @header = @columns;
                } else {
                    my %entry;
                    for (my $column = 0; $column <= $#header; $column++) {
                        $entry{$header[$column]} = $columns[$column];
                    }

                    push(@csv, \%entry);
                }
            }
        }
    }
    close(CSV);

    return @csv;
}

my %args = &getCommandArguments("EVAL_ROOT" => undef, "OUT_ROOT" => undef, "POP" => undef);

my %values;

my @header;
my @counts = &getCSV("$args{'EVAL_ROOT'}.Count_Variants.csv");
my @titv = &getCSV("$args{'EVAL_ROOT'}.Ti_slash_Tv_Variant_Evaluator.csv");

my $varout = "$args{'OUT_ROOT'}.kvp";
open(VAROUT, ">$varout");
foreach my $entry (@counts) {
    my %entry = %$entry;

    print VAROUT "$entry{'jexl_expression'}.$entry{'filter_name'}.$entry{'novelty_name'}.nVariantLoci\t$entry{'nVariantLoci'}\n";

    if ($entry{'jexl_expression'} eq 'Intersection') { $values{'Intersection'} = $entry{'nVariantLoci'}; }
    if ($entry{'jexl_expression'} eq 'UG') { $values{'UG-only'} = $entry{'nVariantLoci'}; }
    if ($entry{'jexl_expression'} eq 'QCALL') { $values{'QCALL-only'} = $entry{'nVariantLoci'}; }
}

foreach my $entry (@titv) {
    my %entry = %$entry;

    print VAROUT "$entry{'jexl_expression'}.$entry{'filter_name'}.$entry{'novelty_name'}.ti/tv_ratio\t$entry{'ti/tv ratio'}\n";
}

my $outpng = "$args{'OUT_ROOT'}.png";
print VAROUT "venndiagram\t$args{'OUT_ROOT'}\n";
print VAROUT "population\t$args{'POP'}\n";
print VAROUT "timestamp\t" . gmtime() . "\n";
close(VAROUT);

my $ug = $values{'UG-only'} + $values{'Intersection'};
my $qcall = $values{'QCALL-only'} + $values{'Intersection'};

my $scaledug = ($ug > $qcall) ? 100 : int($ug*100/$qcall);
my $scaledqcall = ($qcall > $ug) ? 100 : int($qcall*100/$ug);
my $scaledint = ($ug > $qcall) ? int($values{'Intersection'}*100/$ug) : int($values{'Intersection'}*100/$qcall);

my $dm = new DistributedMake("dryRun" => 0);

my $vennCmd = "wget -O $outpng \"http://chart.apis.google.com/chart?cht=v&chd=t:$scaledqcall,$scaledug,0,$scaledint&chco=FF3355,00BBFF&chs=600x480&chdl=QCALL|UnifiedGenotyper&chtt=Callset+concordance+($args{'POP'})&chts=000000,20&chma=100,0,0,0\"";
$dm->addRule($outpng, $varout, $vennCmd);

my $render = "$args{'OUT_ROOT'}.tex";
my $renderCmd = "$ENV{'STING_DIR'}/perl/reports/1KG/RenderTemplate.pl KVP=$varout TEMPLATE=$ENV{'STING_DIR'}/perl/reports/1KG/venn_template.tex TEMPLATE_OUT=$render";
$dm->addRule($render, $outpng, $renderCmd);

my $pdf = "$args{'OUT_ROOT'}.pdf";
my $pdfdir = dirname($args{'OUT_ROOT'});
my $pdfCmd = "pdflatex -output-directory $pdfdir $render";
$dm->addRule($pdf, $render, $pdfCmd);

$dm->execute();
