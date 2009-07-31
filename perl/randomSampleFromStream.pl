#!/usr/bin/perl -w

use strict;
use Getopt::Long;

sub usage {
    print "\nUsage: randomSampleFromStream.pl [-N size] [file1 file2 ...]\n\n";
    print "   Selects a random sample of 'size' elements, without replacement,\n";
    print "   from dataset that represents a union of all input steams. If no\n";
    print "   input streams are specified, then reads STDIN. This script implements\n";
    print "   the standard reservoir sampling algorithm, i.e. it does not preload\n";
    print "   the data into memory and performs selection in one pass.\n\n";
    print "      -N size : optional (default=1), size of the random sample to select.\n\n";

    exit(0);
}



my @selectedLines; # the line we are going to print at the end
my $sampleSize = 1;

my @streams;
my $curr_stream;

sub nextLine {
    my $line = <$curr_stream>;
    
    return $line if ( $line ) ;

    if ( $curr_stream ne "STDIN" && scalar @streams > 0 ) {
        # we are done with the current stream: try opening next one
        until ( $line ) { 
            close $curr_stream;
            if ( scalar @streams > 0 ) {
                my $fname = shift @streams;
                open($curr_stream, "< $fname") or 
                    die("Can not open input file $fname");    
                $line = <$curr_stream>;
            } else {
                last; # no more streams left
            }
        } 
    }
    return $line;
}

my $help = 0;
GetOptions( "N:s" => \$sampleSize,
            "h" => \$help ) or usage();

usage() if ( $help ) ;

if ( scalar(@ARGV) == 0 ) {
    $curr_stream = "STDIN";
} else {
    my $fname = shift @ARGV;
    open($curr_stream, "< $fname") or 
            die("Can not open input file $fname");    
    push @streams, @ARGV;
}


my $line;

for ( my $i = 0 ; $i < $sampleSize; $i++ ) {
    $line = nextLine();
    if ( $line ) {
        push @selectedLines, $line;
    } else {
        # no more lines in the input stream(s)! we got less than sampleSize so far!
        $sampleSize = $i ; # reset sampleSize to the actual number of lines available
        last; 
    }
}


$line = nextLine() if ( $line ) ; # if no more lines left, do not attempt to read
my $index = 0; # where to insert line if selected

my $counter = $sampleSize; # total number of lines read

while ( $line ) {
    $counter++;

    my $prob = $sampleSize/$counter; 

    if ( rand() <= $prob ) {
        # line gets selected
        $index = int ( rand ( $sampleSize ) ) if ( $sampleSize > 1 ); # choose where to insert
        $selectedLines[$index] = $line; # replace old value with newly selected line
    } 
    $line = nextLine();
}


for ( my $i = 0 ; $i < $sampleSize ; $i++ ) { print $selectedLines[$i]; }

