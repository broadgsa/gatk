#!/usr/bin/perl -w

use Getopt::Long;

sub usage {
    print "Usage: perl syncFilesInList.pl\n\t-files <file containing files to sync>\n\t-protocol <protocol to use> [defaults to 'aspera'; can also use 'wget']\n\t[-dry]\n";
    exit(1);
}

my $files = undef;
my $dry;
my $protocol = "aspera";

GetOptions( "files=s"    => \$files,
	    "dry!"       => \$dry,
	    "protocol=s" => \$protocol);

usage() if ( !$files );

open(LIST, "< $files") or die "can't open $files: $!";
while ( <LIST> ) {
    chomp($_);
    if ( $protocol eq "aspera" ) {
	$_ =~ m/data\/(.*)\/alignment.*/;
	$cmd = "./runAspera.pl -source anonftp\@ftp-private.ncbi.nih.gov:/1000genomes/ftp/$_ -dest /humgen/1kg/DCC/ftp/data/$1/alignment/";
	execute($cmd, $dry);
    } elsif ( $protocol eq "wget" ) {
	$cmd = "./runWget.pl -file $_";
	execute($cmd, $dry);
    } else {
	usage();
    }
}
close(LIST);

sub execute {

    my $cmd = $_[0];
    my $dry = $_[1];

    if ($dry) {
	print "$cmd\n";
    } else {
	system($cmd);
    }
}
