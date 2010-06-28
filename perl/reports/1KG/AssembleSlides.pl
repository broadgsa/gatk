#!/usr/bin/perl -w

use strict;
use FindBin;
use lib "$FindBin::Bin/../../";

use StingArgs;

my %args = &getCommandArguments("DIR" => undef, "OUT" => undef);

my $outfile = $args{'OUT'};
chomp(my @files = qx(find $args{'DIR'} -path \*slide\*.pdf | sort -n));

print join("\n", @files) . "\n";

my $cmd = "gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dAutoRotatePages=/All -sOutputFile=$outfile " . join(' ', @files);

system($cmd);
