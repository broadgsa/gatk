#!/usr/bin/perl -w

use strict;
use Data::Dumper;

sub usage {
    print "Usage: $0 <input.vcf> <input.maf> <output.vcf>\n";
    print "This program takes an annotated MAF file and propagates the annotations to the VCF file.\n";
    exit(1);
}

my %args;
($args{'VCF_IN'}, $args{'MAF_IN'}) = @ARGV;

if (!defined($args{'VCF_IN'}) || !defined($args{'MAF_IN'})) {
    &usage();
}

($args{'VCF_OUT'} = $args{'VCF_IN'}) =~ s/\.vcf/.maf_annotated.vcf/;

my %ignoreEntries = (
    'normal_barcode' => 1,
    'tumor_barcode' => 1,
    'build' => 1,
    'tum_allele1' => 1,
    'ref_allele' => 1,
    'tum_allele2' => 1,
    'start' => 1,
    'end' => 1,
);

my %maf;
open(MAF_IN, $args{'MAF_IN'});

chomp(my $mafheader = <MAF_IN>);
chomp($mafheader);
my @mafheader = split(/\s+/, $mafheader);

while (my $mafline = <MAF_IN>) {
    chomp($mafline);
    my @mafline = split(/\s+/, $mafline);

    my %mafentry;
    for (my $i = 0; $i <= $#mafheader; $i++) {
        $mafentry{$mafheader[$i]} = $mafline[$i];
    }

    my $locus = "$mafentry{'chr'}:$mafentry{'start'}";
    if (exists($mafentry{$locus})) {
        print "Locus $locus already spoken for!\n";
        exit(1);
    }

    $maf{$locus} = \%mafentry;
}
close(MAF_IN);

open(VCF_OUT, ">$args{'VCF_OUT'}");
open(VCF_IN, $args{'VCF_IN'});

my @vcfheader;
while (my $vcfline = <VCF_IN>) {
    chomp($vcfline);

    if ($vcfline =~ /##/) {
        print VCF_OUT "$vcfline\n";
    } elsif ($vcfline =~ /#CHROM/) {

        print VCF_OUT "##source=AnnotateVCFwithMAF\n";
        print VCF_OUT "##INFO=<ID=cDNAchange,Number=1,Type=String,Description=\"cDNAchange\">\n";
        print VCF_OUT "##INFO=<ID=classification,Number=1,Type=String,Description=\"classification\">\n";
        print VCF_OUT "##INFO=<ID=codonchange,Number=1,Type=String,Description=\"codonchange\">\n";
        print VCF_OUT "##INFO=<ID=gene,Number=1,Type=String,Description=\"gene\">\n";
        print VCF_OUT "##INFO=<ID=genomechange,Number=1,Type=String,Description=\"genomechange\">\n";
        print VCF_OUT "##INFO=<ID=proteinchange,Number=1,Type=String,Description=\"proteinchange\">\n";
        print VCF_OUT "##INFO=<ID=strand,Number=1,Type=String,Description=\"strand\">\n";
        print VCF_OUT "##INFO=<ID=transcript,Number=1,Type=String,Description=\"transcript\">\n";
        print VCF_OUT "##INFO=<ID=type,Number=1,Type=String,Description=\"type\">\n";
        print VCF_OUT "$vcfline\n";

        $vcfline =~ s/#//g;
        @vcfheader = split(/\s+/, $vcfline);
    } else {
        my @vcfentry = split(/\s+/, $vcfline);
        my %vcfentry;

        for (my $i = 0; $i <= $#vcfheader; $i++) {
            $vcfentry{$vcfheader[$i]} = $vcfentry[$i];
        }

        my $locus = "$vcfentry{'CHROM'}:$vcfentry{'POS'}";
        if (!exists($maf{$locus})) {
            for (my $i = 1; $i <= 10; $i++) {
                $locus = "$vcfentry{'CHROM'}:" . ($vcfentry{'POS'}-$i);

                last if (exists($maf{$locus}));
            }
        }
        my %mafentry = %{$maf{$locus}};

        my @info = split(/;/, $vcfentry{'INFO'});
        my %info;
        foreach my $info (@info) {
            my ($key, $value) = split(/=/, $info);

            $info{$key} = $value;
        }

        foreach my $mafkey (sort { $a cmp $b } keys(%mafentry)) {
            if (!$ignoreEntries{$mafkey}) {
                $info{$mafkey} = $mafentry{$mafkey};
            }
        }

        my @newinfo;
        foreach my $infokey (sort { $a cmp $b } keys(%info)) {
            if (!defined($info{$infokey})) {
                #print "$infokey is missing\n";
            } else {
                push(@newinfo, "$infokey=$info{$infokey}");
            }
        }
        $vcfentry{'INFO'} = join(";", @newinfo);

        my @newvcfline;
        for (my $i = 0; $i <= $#vcfheader; $i++) {
            push(@newvcfline, $vcfentry{$vcfheader[$i]});
        }

        print VCF_OUT join("\t", @newvcfline) . "\n";
    }
}

close(VCF_IN);
close(VCF_OUT);
