#!/usr/bin/perl -w

# Runs Wget to pull down a file

use strict;
use Getopt::Long;

my $file = undef;
GetOptions( "file=s" => \$file);

if ( !$file) {
    print "Usage: runWget.pl\n\t-file \t<file>\n";
    exit(1);
}

chomp($file);
my $cmd = "wget -O /humgen/1kg/DCC/ftp/$file ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/$file";
print "$cmd\n";
system($cmd);
