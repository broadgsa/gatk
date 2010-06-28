#!/usr/bin/perl -w

use strict;
use lib "$ENV{'STING_DIR'}/perl";
use StingArgs;

my %args = &getCommandArguments("VCF" => undef, "OUT" => "/dev/stdout");

my %table;
my %sets;

open(VCF, $args{'VCF'});

my @header;
while (my $line = <VCF>) {
    chomp($line);
    if ($line =~ /#/) {
        if ($line =~ /CHROM/) {
            $line =~ s/#//g;
            @header = split(/\s+/, $line);
        }
    } else {
        my @columns = split(/\s+/, $line);

        if ($columns[0] eq '1') {
            my %entry;
            for (my $i = 0; $i <= $#columns; $i++) {
                $entry{$header[$i]} = $columns[$i];
            }

            my ($set) = $entry{'INFO'} =~ /set=(\w+)/;
            $sets{$set} = 1;

            for (my $i = 9; $i <= $#header; $i++) {
                if ($entry{$header[$i]} =~ /0[\\\|\/]1/ ||
                    $entry{$header[$i]} =~ /1[\\\|\/]0/ ||
                    $entry{$header[$i]} =~ /1[\\\|\/]1/) {

                    ${$table{$header[$i]}}{$set}++;
                    $sets{$set} = 1;
                }
            }
        }
    }
}
close(VCF);

open(OUT, ">$args{'OUT'}");

print OUT "sample\t" . join("\t", sort { $a cmp $b } keys(%sets)) . "\n";

foreach my $sample (sort { $a cmp $b } keys(%table)) {
    print OUT "$sample";

    foreach my $set (sort { $a cmp $b } keys(%sets)) {
        my $value = 0;
        if (exists(${$table{$sample}}{$set})) {
            $value = ${$table{$sample}}{$set};
        }
        
        print OUT "\t$value";
    }

    print OUT "\n";
}

close(OUT);
