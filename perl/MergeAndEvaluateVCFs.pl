#!/usr/bin/perl -w

use strict;
use Getopt::Long;

sub getArgs {
    my %doc = (
        'GATK'          => "Path to GenomeAnalysisTK.jar (usually in your Sting/dist/ directory)",
        'VARIANTREPORT' => "Path to VariantReport.R (usually in your Sting/R/VariantReport/ directory)",
        'REFERENCE'     => "Path to human reference sequence",
        'DBSNP'         => "Path to dbsnp rod",
        'VCF1'          => "Path to first VCF file",
        'VCF1NAME'      => "Name of first VCF file",
        'VCF2'          => "Path to second VCF file",
        'VCF2NAME'      => "Name of second VCF file",
        'VCFOUT'        => "Output root for resultant merged VCF, the eval file, and PDF",
    );

    my %args = (
        'GATK'          => undef,
        'VARIANTREPORT' => undef,
        'REFERENCE'     => undef,
        'DBSNP'         => undef,
        'VCF1'          => undef,
        'VCF1NAME'      => undef,
        'VCF2'          => undef,
        'VCF2NAME'      => undef,
        'VCFOUT'        => undef,
    );

    my $result = GetOptions(
        "GATK=s"          => \$args{'GATK'},
        "VARIANTREPORT=s" => \$args{'VARIANTREPORT'},
        "REFERENCE=s"     => \$args{'REFERENCE'},
        "DBSNP=s"         => \$args{'DBSNP'},
        "VCF1=s"          => \$args{'VCF1'},
        "VCF1NAME=s"      => \$args{'VCF1NAME'},
        "VCF2=s"          => \$args{'VCF2'},
        "VCF2NAME=s"      => \$args{'VCF2NAME'},
        "VCFOUT=s"        => \$args{'VCFOUT'},
    );

    my @undefinedArgs;
    foreach my $key (sort { $a cmp $b } keys(%args)) {
        if (!defined($args{$key})) {
            push(@undefinedArgs, $key);
        }
    }

    if (scalar(@undefinedArgs) > 0) {
        print "Error: there were some undefined arguments:\n";

        foreach my $undefinedArg (@undefinedArgs) {
            print "\t-$undefinedArg : $doc{$undefinedArg}\n";
        }

        exit(-1);
    }

    return %args;
}

sub runCommand {
    my ($output, $cmd) = @_;

    #if (!-e $output) {
        system($cmd);
    #}
}

my %args = &getArgs();

my $combineVariantsCmd = "java -jar $args{'GATK'} -T CombineVariants -R $args{'REFERENCE'} -B:$args{'VCF1NAME'},VCF $args{'VCF1'} -B:$args{'VCF2NAME'},VCF $args{'VCF2'} -variantMergeOptions UNION -priority $args{'VCF1NAME'},$args{'VCF2NAME'} -o $args{'VCFOUT'} -l INFO";
&runCommand($args{'VCFOUT'}, $combineVariantsCmd);

my $variantEvalCmd = "java -jar $args{'GATK'} -T VariantEval -R $args{'REFERENCE'} -D $args{'DBSNP'} -B:eval,VCF $args{'VCFOUT'} -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"$args{'VCF1NAME'}\"' -selectName $args{'VCF1NAME'} -select 'set==\"filterIn$args{'VCF1NAME'}-$args{'VCF2NAME'}\"' -selectName In$args{'VCF2NAME'}-FilteredIn$args{'VCF1NAME'} -select 'set==\"$args{'VCF2NAME'}\"' -selectName $args{'VCF2NAME'} -select 'set==\"filterIn$args{'VCF2NAME'}-$args{'VCF1NAME'}\"' -selectName In$args{'VCF1NAME'}-FilteredIn$args{'VCF2NAME'} -select 'set==\"FilteredInAll\"' -selectName FilteredInAll -reportType R -reportLocation $args{'VCFOUT'}.eval";
&runCommand("$args{'VCFOUT'}.eval", $variantEvalCmd);

my $variantReportCmd = "Rscript $args{'VARIANTREPORT'} -title 'Automated Variant Report' -author '$ENV{'USER'}' -evalRoot $args{'VCFOUT'}.eval -plotOut $args{'VCFOUT'}.eval.pdf";
&runCommand("$args{'VCFOUT'}.eval.pdf", $variantReportCmd);
