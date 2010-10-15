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
my $cacheSize = 50000;
my $tmp = "/tmp";
GetOptions( "vcf=s" => \$in,
	    "gatk=s" => \$gatk,
	    "chain=s" => \$chain,
	    "newRef=s" => \$newRef,
	    "oldRef=s" => \$oldRef,
            "out=s" => \$out,
	    "cache=s" => \$cacheSize,
	    "tmp=s" => \$tmp);

if ( !$in || !$gatk || !$chain || !$newRef || !$oldRef || !$out ) {
    print "Usage: liftOverVCF.pl\n\t-vcf \t\t<input vcf>\n\t-gatk \t\t<path to gatk trunk>\n\t-chain \t\t<chain file>\n\t-newRef \t<path to new reference prefix; we will need newRef.dict, .fasta, and .fasta.fai>\n\t-oldRef \t<path to old reference prefix; we will need oldRef.fasta>\n\t-out \t\t<output vcf>\n\t-cache \t\t<cache size used for on-the-fly sorting; defaults to 50000>\n\t-tmp \t\t<temp file location; defaults to /tmp>\n";
    print "Example: ./liftOverVCF.pl\n\t-vcf /humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/1kg_snp_validation/all_validation_batches.b36.vcf\n\t-chain b36ToHg19.broad.over.chain\n\t-out lifted.hg19.vcf\n\t-gatk /humgen/gsa-scr1/ebanks/Sting_dev\n\t-newRef /seq/references/Homo_sapiens_assembly19/v0/Homo_sapiens_assembly19\n\t-oldRef /humgen/1kg/reference/human_b36_both\n\t-cache 60000\n";
    exit(1);
}

# generate a random number
my $random_number = rand();
my $tmp_prefix = "$tmp/$random_number";
print "Writing temporary files to prefix: $tmp_prefix\n";
my $lifted_vcf = "$tmp_prefix.lifted.vcf";

# lift over the file
print "Lifting over the vcf...";
my $cmd = "java -jar $gatk/dist/GenomeAnalysisTK.jar -T LiftoverVariants -R $oldRef.fasta -B:variant,vcf $in -o $lifted_vcf -chain $chain -dict $newRef.dict -cache $cacheSize";
system($cmd);

# Filter the VCF for bad records
print "\nFixing/removing bad records...\n";
$cmd = "java -jar $gatk/dist/GenomeAnalysisTK.jar -T FilterLiftedVariants -R $newRef.fasta -B:variant,vcf $lifted_vcf -o $out";
system($cmd);

# clean up
unlink $lifted_vcf;

print "\nDone!\n";
