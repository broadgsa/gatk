#!/usr/bin/perl -w

# Runs the liftover tool on a VCF and properly handles the output

use strict;
use Getopt::Long;

my $in = undef;
my $gatk = undef;
my $chain = undef;
my $newRef = undef;
my $oldRef = undef;
my $out = undef;
my $tmp = "/tmp";
my $recordOriginalLocation = 0;
GetOptions( "vcf=s" => \$in,
	    "gatk=s" => \$gatk,
	    "chain=s" => \$chain,
	    "newRef=s" => \$newRef,
	    "oldRef=s" => \$oldRef,
            "out=s" => \$out,
	    "tmp=s" => \$tmp,
	    "recordOriginalLocation" => \$recordOriginalLocation);

if ( !$in || !$gatk || !$chain || !$newRef || !$oldRef || !$out ) {
    print "Usage: liftOverVCF.pl\n\t-vcf \t\t<input vcf>\n\t-gatk \t\t<path to gatk trunk>\n\t-chain \t\t<chain file>\n\t-newRef \t<path to new reference prefix; we will need newRef.dict, .fasta, and .fasta.fai>\n\t-oldRef \t<path to old reference prefix; we will need oldRef.fasta>\n\t-out \t\t<output vcf>\n\t-tmp \t\t<temp file location; defaults to /tmp>\n\t-recordOriginalLocation \t\t<Should we record what the original location was in the INFO field?; defaults to false>\n";
    print "Example: ./liftOverVCF.pl\n\t-vcf /humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/1kg_snp_validation/all_validation_batches.b36.vcf\n\t-chain b36ToHg19.broad.over.chain\n\t-out lifted.hg19.vcf\n\t-gatk /humgen/gsa-scr1/ebanks/Sting_dev\n\t-newRef /seq/references/Homo_sapiens_assembly19/v0/Homo_sapiens_assembly19\n\t-oldRef /humgen/1kg/reference/human_b36_both\n";
    exit(1);
}

# generate a random number
my $random_number = rand();
my $tmp_prefix = "$tmp/$random_number";
print "Writing temporary files to prefix: $tmp_prefix\n";
my $unsorted_vcf = "$tmp_prefix.unsorted.vcf";

# lift over the file
print "Lifting over the vcf...";
my $cmd = "java -jar $gatk/dist/GenomeAnalysisTK.jar -T LiftoverVariants -R $oldRef.fasta -V:variant $in -o $unsorted_vcf -chain $chain -dict $newRef.dict";
if ($recordOriginalLocation) {
  $cmd .= " -recordOriginalLocation";
}
system($cmd) == 0 or quit("The liftover step failed.  Please correct the necessary errors before retrying.");

# we need to sort the lifted over file now
print "\nRe-sorting the vcf...\n";
my $sorted_vcf = "$tmp_prefix.sorted.vcf";
open(SORTED, ">$sorted_vcf") or die "can't open $sorted_vcf: $!";

# write the header
open(UNSORTED, "< $unsorted_vcf") or die "can't open $unsorted_vcf: $!";
my $inHeader = 1;
while ( $inHeader == 1 ) {
    my $line = <UNSORTED>;
    if ( $line !~ m/^#/ ) {
	$inHeader = 0;
    } else {
	print SORTED "$line";
    }
}
close(UNSORTED);
close(SORTED);

$cmd = "grep \"^#\" -v $unsorted_vcf | sort -n -k2 -T $tmp | $gatk/public/perl/sortByRef.pl --tmp $tmp - $newRef.fasta.fai >> $sorted_vcf";
system($cmd) == 0 or quit("The sorting step failed.  Please correct the necessary errors before retrying.");

# Filter the VCF for bad records
print "\nFixing/removing bad records...\n";
$cmd = "java -jar $gatk/dist/GenomeAnalysisTK.jar -T FilterLiftedVariants -R $newRef.fasta -V:variant $sorted_vcf -o $out";
system($cmd) == 0 or quit("The filtering step failed.  Please correct the necessary errors before retrying.");

# clean up
unlink $unsorted_vcf;
unlink $sorted_vcf;
my $sorted_index = "$sorted_vcf.idx";
unlink $sorted_index;

print "\nDone!\n";

sub quit {
    print "\n$_[0]\n";
    exit(1);
}
