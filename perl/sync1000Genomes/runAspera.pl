#!/usr/bin/perl -w

# Runs Aspera to pull down files

use strict;
use Getopt::Long;

my $source = undef;
my $dest = ".";
GetOptions( "source=s" => \$source,
	    "dest=s" => \$dest);

if ( !$source) {
    print "Usage: runAspera.pl\n\t-source \t<ftp source>\n\t-dest \t\t<local destination; defaults to '.'>\n";
    exit(1);
}

my $cmd = "ascp -i /opt/aspera/etc/asperaweb_id_dsa.putty -k2 -QTr -l2G -d -v $source $dest";
system($cmd);
