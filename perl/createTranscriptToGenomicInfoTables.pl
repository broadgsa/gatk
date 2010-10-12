#!/usr/bin/perl -w

# Runs the TranscriptToGenomicInfo tool to create big tables

use strict;
use Getopt::Long;

my $in = undef;
my $gatk = undef;
my $ref = undef;
my $out = undef;
my $tmp = "/tmp";
GetOptions( "transcript=s" => \$in,
	    "gatk=s" => \$gatk,
	    "ref=s" => \$ref,
            "out=s" => \$out,
	    "tmp=s" => \$tmp);

if ( !$in || !$gatk || !$ref || !$out ) {
    print "Usage: createTranscriptToGenomicInfoTables.pl\n\t-transcript \t<raw transcript input table>\n\t-gatk \t\t<path to gatk trunk>\n\t-ref \t\t<path to reference.fasta>\n\t-out \t\t<genomic output table>\n\t-tmp \t\t<temp file location; needs potentially lots of space; defaults to /tmp>\n";
    print "Example: ./createTranscriptToGenomicInfoTables.pl\n\t-transcript test.foo\n\t-ref /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta\n\t-gatk /humgen/gsa-scr1/ebanks/Sting_dev\n\t-out test.bar\n\t-tmp /broad/shptmp\n";
    exit(1);
}

# generate a random number
my $random_number = rand();
my $tmp_prefix = "$tmp/$random_number";
print "Writing temporary files to prefix: $tmp_prefix\n";
my $unsorted_table = "$tmp_prefix.unsorted.table";

# convert the file
print "Converting the transcript table...";
my $cmd = "java -jar $gatk/dist/GenomeAnalysisTK.jar -T TranscriptToGenomicInfo -R $ref -B:transcripts,AnnotatorInputTable $in -o $unsorted_table -n name,proteinID";
system($cmd);

# we need to sort the converted table now
print "\nRe-sorting the table...\n";
open(SORTED, ">$out") or die "can't open $out: $!";

# write the header
open(UNSORTED, "< $unsorted_table") or die "can't open $unsorted_table: $!";
my $line = <UNSORTED>;
print SORTED "$line";
close(UNSORTED);

$cmd = "grep haplotypeReference -v $unsorted_table | sort -n -k2 -T $tmp | $gatk/perl/sortByRef.pl --tmp $tmp - $ref.fai";
print SORTED `$cmd`;
close(SORTED);

# clean up
unlink $unsorted_table;

print "\nDone!\n";
