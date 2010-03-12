package DataTable;

use strict;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
require Exporter;
 
@ISA = qw(Exporter AutoLoader);

@EXPORT = qw( readTable );

use Data::Dumper;

sub readTable {
    my ($file, %overrides) = @_;
    my %options = (
        'file' => $file,
        'key' => undef,
        'header' => 0,
        'delimiter' => '\s+',
        'filter' => '#',
        %overrides,
    );

    open(FILE, $options{'file'});

    my @header = undef;
    if ($options{'header'} == 1) {
        chomp(my $headerline = <FILE>);
        $headerline =~ s/#//g;
        @header = split(/$options{'delimiter'}/, $headerline);
    } elsif (ref($options{'header'}) eq 'ARRAY') {
        @header = @{$options{'header'}};
    }

    chomp(my @lines = grep { $_ !~ /$options{'filter'}/ } <FILE>);

    my @table;
    my %table;
    for (my $lineIndex = 0; $lineIndex <= $#lines; $lineIndex++) {
        my $line = $lines[$lineIndex];

        my %fieldHash;
        my @fields = split(/$options{'delimiter'}/, $line);

        for (my $fieldIndex = 0; $fieldIndex <= $#fields; $fieldIndex++) {
            $fieldHash{$header[$fieldIndex]} = $fields[$fieldIndex];
        }

        $fieldHash{'_linenum'} = $lineIndex;
        $fieldHash{'_line'} = $line;

        if (!defined($options{'key'})) {
            push(@table, \%fieldHash);
        } else {
            my $key = ($options{'key'} =~ /^\d+/) ? $fields[$options{'key'}] : $fieldHash{$options{'key'}};

            push(@{$table{$key}}, \%fieldHash);
        }
    }

    return (!defined($options{'key'})) ? @table : %table;
}

1;
