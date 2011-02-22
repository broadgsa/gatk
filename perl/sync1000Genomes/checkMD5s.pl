#!/usr/bin/perl -w

use Getopt::Long;

sub usage {
    print "Usage: perl checkMD5s.pl\n\t-ai <alignment.index to check>\n\t-o <file to store results>\n";
    exit(1);
}


my $ai = undef;
my $out = undef;
GetOptions( "ai=s" => \$ai,
	    "o=s" => \$out);

usage() if ( !$ai || !$out );

open(OUT, "> $out") or die "can't open $out: $!";

open(LIST, "< $ai") or die "can't open $ai: $!";
while ( <LIST> ) {
    @pieces = split(' ', $_);
    if ( @pieces == 6 ) {
	check($pieces[0], $pieces[1]);
	check($pieces[2], $pieces[3]);
	check($pieces[4], $pieces[5]);
    }
}

close(LIST);
close(OUT);

sub check {

    my $file = $_[0];
    my $target = $_[1];

    print "Checking /humgen/1kg/DCC/ftp/$file\n";
    @md5 = split(' ', `md5sum /humgen/1kg/DCC/ftp/$file`);
    if ( $md5[0] ne $target ) {
	print OUT "$file\t$md5[0]\t$target\n";
    }
}
