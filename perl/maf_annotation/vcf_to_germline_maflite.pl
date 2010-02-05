#!/util/bin/perl
$|++;
use strict;

if (scalar(@ARGV) != 1 ) {
    die("usage vcf_to_germline_maflite.pl <normalAlias>\n");
}
 
my ($normalAlias) = @ARGV;
my $line;
while (defined($line = <STDIN>)) {
    if ($line =~ /#/) { next; }

    my ($chrom, $pos, $rsName, $ref, $alt, $qual, $filter, $info, $format, $sampleA) = split("\t", $line);

    my ($genoNum) = split(":", $sampleA);
    my ($geno1, $geno2) = split("/",$genoNum);
       
    my $allele1 = "N";
    my $allele2 = "N";
    if ( $geno1 eq "0") { 
        $allele1 = $ref;
    } elsif ($geno1 eq "1") {
        $allele1 = $alt;
    } else {
        die("only handled single alt currently:\n$line \n");
    }

    if ( $geno2 eq "0") { 
        $allele2 = $ref;
    } elsif ($geno2 eq "1") {
        $allele2 = $alt;
    } else {
        die("only handled single alt currently:\n$line\n");
    }

    my $ncbiBuild = "36";
    print join("\t", ($ncbiBuild, $chrom, $pos, $pos, $ref, $allele1, $allele2, $normalAlias, $normalAlias)) . "\n";
}
