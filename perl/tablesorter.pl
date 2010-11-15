#!/usr/bin/perl -w

# sorts by reference

use strict;
use Getopt::Long;

my $in = undef;
my $gatk = undef;
my $ref = undef;
my $out = undef;
my $tmp = "/tmp";
my $contig = undef;
my $startpos =undef;
my $mark =undef;

GetOptions( "in=s" => \$in,
            "gatk=s" => \$gatk,
            "ref=s" => \$ref,
            "out=s" => \$out,
            "tmp=s" => \$tmp,
            "startpos=s" => \$startpos,
            "contig=s" => \$contig,
            "headermarker=s" => \$mark);

	    if ( !$in || !$gatk || !$ref || !$out ) {
		print "Usage: tablesorter.pl\n\t-in \t\t\t<unsorted input table>\n\t-gatk \t\t\t<path to gatk trunk>\n\t-ref \t\t\t<path to reference.fasta>\n\t-out \t\t\t<sorted output table>\n\t-tmp \t\t\t<temp file location; needs potentially lots of space; defaults to /tmp>\n\t-startpos \t\t<Which column is your start position listed in?>\n\t-contig \t\t<Which column is your contig listed in?>\n\t-headermarker \t\t<Give a field name unique to your header (not found in your table).>\n";
		print "Example: ./tablesorter.pl\n\t-in test.foo\n\t-ref /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta\n\t-gatk /humgen/gsa-s r1/ebanks/Sting_dev\n\t-out test.bar\n\t-tmp /broad/shptmp\n\t-startpos /2\n\t-contig 1\n\t-headermarker haplotype \n";
		exit(1);
	    }


# we need to sort the converted table now
	    print "\nSorting the table...\n";
	    open(SORTED, ">$out") or die "can't open $out: $!";

# write the header
	    open(UNSORTED, "< $in") or die "can't open $in: $!";
	    my $line = <UNSORTED>;
	    print SORTED "$line";
	    close(UNSORTED);
	    close(SORTED);

print "$line is the header";
	    my $cmd = "grep $mark -v $in | sort -n -k$startpos -T $tmp | $gatk/perl/sortByRef.pl --k $contig --tmp $tmp - $ref.fai >> $out";
	    system($cmd);
# clean upunlink $unsorted_table;

	    print "\nDone!\n";
