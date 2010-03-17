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

    my ($chrom, $pos, $rsName, $ref, $alt) = split("\t", $line);

    my $ncbiBuild = "36";
    print join("\t", ($ncbiBuild, $chrom, $pos, $pos, $ref, $alt, $alt, $normalAlias, $normalAlias)) . "\n";
}
