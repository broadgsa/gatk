#!/usr/bin/perl -w

use strict;
use lib "$ENV{'STING_DIR'}/perl";
use StingArgs;
use Data::Dumper;

my %args = &getCommandArguments("KVP" => undef, "TEMPLATE" => undef, "TEMPLATE_OUT" => undef);

my %vars;
open(KVP, $args{'KVP'});
while (my $kvpline = <KVP>) {
    chomp($kvpline);

    my ($key, $value) = split(/\t+/, $kvpline);
    $vars{$key} = $value;
}
close(KVP);

open(TEMPLATE, $args{'TEMPLATE'});
open(TEMPLATE_OUT, ">$args{'TEMPLATE_OUT'}");
while (my $line = <TEMPLATE>) {
    chomp($line);

    while ($line =~ /\$(.+?)\$/) {
        if (exists($vars{$1})) {
            my ($key, $value) = ($1, $vars{$1});
            $line =~ s/\$$key\$/$value/;
        } else {
            $line =~ s/\$$1\$/unknown/;
        }
    }

    print TEMPLATE_OUT "$line\n";
}
close(TEMPLATE_OUT);
close(TEMPLATE);
