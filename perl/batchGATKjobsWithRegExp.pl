#!/usr/bin/perl -w

# Runs a given (GATK) command as specified, except that a single '*'
# in the -I argument field is replaced with the appropriate shell
# expansion and all other instances of the '*' in other arguments are
# replaced with the correct expansion accordingly.
# IMPORTANT: remember to surround your command with quotes (\"\") so
# that the shell doesn't try to expand out your regular expression!

use strict;
use Getopt::Long;

my $dry;
my $inArg = "I";
GetOptions( "dry!" => \$dry,
            "inArg=s" => \$inArg);

if (scalar(@ARGV) != 1) {
    print "Usage: batchGATKjobsWithRegExp\n\t[-inArg <the GATK argument flag for input; default=I>]\n\t[-dry]\n\t\"GATK command\"\n";
    exit(1);
}

my @args = split(/ /, $ARGV[0]);
chomp(@args);

my $argcount = scalar(@args);
my $IargIdx = undef;
for (my $i = 0; $i < $argcount; $i++) {
    if ($args[$i] eq "-$inArg") {
	$IargIdx = $i + 1;
    }
}
my $Iarg = "$args[$IargIdx]\n";
if ($Iarg =~ m/(.*),(.*),(.*)/) {
    $Iarg = $3;
}
my @matches = glob($Iarg);
$Iarg =~ s/\*/(.*)/;
chomp($Iarg);

foreach my $match (@matches) {
    $match =~ m/$Iarg/;
    my $piece = $1;

    my $cmd = "";
    for (my $i = 0; $i < $argcount; $i++) {
	my $arg = $args[$i];
	$arg =~ s/\*/$piece/;
	$cmd .= "$arg ";
    }

    if ($dry) {
	print "$cmd\n";
    } else {
	system($cmd);
    }
}
