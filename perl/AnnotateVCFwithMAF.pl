#!/usr/bin/perl -w

use strict;
use lib "$ENV{'STING_DIR'}/perl";
use StingArgs;
use Data::Dumper;

my %args = &getCommandArguments("VCF_IN" => undef, "MAF_IN" => undef, "VCF_OUT" => "/dev/stdout");

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
        print Dumper(\%mafentry);
        exit(1);
    }

    $maf{$locus} = \%mafentry;

    #print Dumper(\%mafentry);
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
        print VCF_OUT "##INFO=cDNAchange,1,String,cDNAchange\n";
        print VCF_OUT "##INFO=classification,1,String,classification\n";
        print VCF_OUT "##INFO=codonchange,1,String,codonchange\n";
        print VCF_OUT "##INFO=gene,1,String,gene\n";
        print VCF_OUT "##INFO=genomechange,1,String,genomechange\n";
        print VCF_OUT "##INFO=proteinchange,1,String,proteinchange\n";
        print VCF_OUT "##INFO=strand,1,String,strand\n";
        print VCF_OUT "##INFO=transcript,1,String,transcript\n";
        print VCF_OUT "##INFO=type,1,String,type\n";
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
            push(@newinfo, "$infokey=$info{$infokey}");
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
