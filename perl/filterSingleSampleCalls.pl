#!/usr/bin/perl

use Getopt::Long;
use strict;


sub usage {

    my $message = shift;

    print "\nERROR:\n$message\n\n";

    print "Usage:\n\n";
    print "  filterSingleSampleCalls --calls FILEIN [--max_cons_av_mm N1] [--max_ref_av_mm N2] [--max_cons_nqs_av_mm N3] [--min_cons_nqs_av_qual N4] [--min_ref_nqs_av_qual N5] [--mode MODE]\n\n";
    print "     FILEIN   File to read and apply filter to. If \"-\" then read from stdin.\n";
    print "         N1   max. average number of mismatches per (consensus) indel-containing read.\n";
    print "              If the number is greater then N1, indel will be discarded/marked.\n";
    print "         N2   max. average number of mismatches per reference-matching read.\n";
    print "              If the number is greater then N2, indel will be discarded/marked.\n";
    print "         N3   max. average mismatch rate in NQS window around the indel, across all indel-containing read.\n";
    print "              If the number is greater then N3, indel will be discarded/marked.\n";
    print "         N4   min. average base quality in all indel supporting reads in the nqs window around the indel.\n";
    print "              If the average base quality is less than N4, the indel will be discarded/marked.\n";
    print "         N5   min. average base quality in all reference supporting reads in the nqs window around the indel.\n";
    print "              If the average base quality is less than N5, the indel will be discarded/marked.\n";
    print "       MODE   If (any prefix of) ANNOTATE, then indel calls not passing any of the filters will be still printed into stdout with\n";
    print "              additional AUTOFILTER_* tags added, specifying what cutoff(s) were not passed.\n";
    print "              If (any prefix of) DISCARD, then only indels passing all the filters will be printed into stdout, the rest\n";
    print "              will be discarded (default).\n\n";

    exit(1);

}

my $calls = "";
my $tolerate_old_calls = 1;
my $cons_av_mm_cutoff = 100000;
my $ref_av_mm_cutoff = 100000;
my $cons_av_nqs_mm_cutoff = 100000;
my $ref_nqs_av_qual_cutoff = 0;
my $cons_nqs_av_qual_cutoff = 0;
my $mode_arg = "";
my $mode = 1; # "discard"

GetOptions("calls:s" => \$calls,
           "max_cons_av_mm:f" => \$cons_av_mm_cutoff,
           "max_ref_av_mm:f" => \$ref_av_mm_cutoff,
           "max_cons_nqs_av_mm:f" => \$cons_av_nqs_mm_cutoff,
           "min_ref_nqs_av_qual:f" => \$ref_nqs_av_qual_cutoff,
           "min_cons_nqs_av_qual:f" => \$cons_nqs_av_qual_cutoff,
           "mode=s" => \$mode_arg
    ) 
    or usage("Can not parse argument string");

usage ("--calls argument must be specified") if ( $calls eq "" ) ;
usage ("--mode argument must be specified (unique prefix of: ANNOTATE or DISCARD)") if ( $mode_arg eq "");

if ( "ANNOTATE" =~ /^$mode_arg/ ) {
    $mode = 0;
} elsif ( "DISCARD" =~ /^$mode_arg/ ) {
    $mode=1;
} else {
    die("Unrecognized value specified in --mode argument");
}

my $input_stream;

if ( $calls eq "-" ) {
    $input_stream = "STDIN";
} else {
    open ( $input_stream, "< $calls") or 
        die("Can not open input file $calls: $!");
}

my $id_counter = 0;

while( <$input_stream> ) {

    chomp;

    $id_counter++;
    my $annotation="";

    next if ( $_ =~ /_TOO_DIRTY/ );

#    print $_,"\n";
#    next;

    my $cons_cnt = $1;
    my $indel_cnt = $2;
    my $cov = $3;

    if ( $_ =~ /\sOBS_COUNTS\[C\/A\/R\]:(\d+)\/(\d+)\/(\d+)\s/ ) {
        $cons_cnt = $1;
        $indel_cnt = $2;
        $cov = $3;
    } else {
        if ( $tolerate_old_calls != 0 ) {
            print "$_\n";
            next;
        } else {
            die("NO OBS_COUNTS in\n$_\n");
        }
    }


    die ("NO AV_MM MATCH in\n$_\n") if ( $_ !~ /\sAV_MM\[C\/R\]:([\+-\d\.]+)\/([\+-\d\.]+)\s/ ) ;

    my $cons_av_mm = $1;
    my $ref_av_mm = $2;

    die("NO AV_MAPQ MATCH in\n$_\n") if ( $_ !~ /\sAV_MAPQ\[C\/R\]:([\+-\d\.]+)\/([\+-\d\.]+)\s/ ) ;

    my $cons_av_mapq = $1;
    my $ref_av_mapq = $2;

    die("NO NQS_MM_RATE MATCH in\n$_\n") if ( $_ !~ /\sNQS_MM_RATE\[C\/R\]:([\d\.]+)\/([\d\.]+)\s/ ) ;

    my $cons_nqs_mm_rate = $1;
    my $ref_nqs_mm_rate = $2;

    die( "NO NQS_AV_QUAL MATCH in\n$_\n") if ( $_ !~ /\sNQS_AV_QUAL\[C\/R\]:([\d\.]+)\/([\d\.]+)\s/ ) ;

    my $cons_nqs_av_qual = $1;
    my $ref_nqs_av_qual = $2;


    if ( $cons_av_mm < 0 ) {
        print STDERR "WARNING: negative mismatch rate in consensus supporting reads:\n";
        print STDERR "$_\n";
        next;
    }
    if ( $ref_av_mm < 0 ) {
        print STDERR "WARNING: negative mismatch rate in reference supporting reads:\n";
        print STDERR "$_\n";
        next;
    }
    if ( $cons_av_mapq < 0 ) {
        print STDERR "WARNING: negative average mapping quality in consensus supporting reads:\n";
        print STDERR "$_\n";
        next;
    }
    if ( $ref_av_mapq < 0 ) {
        print STDERR "WARNING: negative average mapping quality in reference supporting reads:\n";
        print STDERR "$_\n";
        next;
    }

    if ( $cons_av_mm > $cons_av_mm_cutoff ) {
        # filter indel out: alignments for indel-containing reads are too messy
        if ( $mode == 0 ) { # annotate
            $annotation .= "\tAUTOFILTER_CONS_AV_MM_$cons_av_mm_cutoff";
        } else {
            next; # discard
        }
    }

    if ( $ref_av_mm > $ref_av_mm_cutoff ) {
        # filter indel out: alignments for reference-matching reads are too messy
        if ( $mode == 0 ) { # annotate
            $annotation .= "\tAUTOFILTER_REF_AV_MM_$ref_av_mm_cutoff";
        } else {
            next; # discard
        }
    }


    if ( $cons_nqs_av_qual < $cons_nqs_av_qual_cutoff ) {
        # filter indel out: alignments for indel-containing reads are too messy
        if ( $mode == 0 ) { # annotate
            $annotation .= "\tAUTOFILTER_CONS_NQS_AV_QUAL_$cons_nqs_av_qual_cutoff";
        } else {
            next; # discard
        }
    }

    if ( $ref_nqs_av_qual < $ref_nqs_av_qual_cutoff ) {
        # filter indel out: alignments for reference-matching reads are too messy
        if ( $mode == 0 ) { # annotate
            $annotation .= "\tAUTOFILTER_REF_NQS_AV_QUAL_$ref_nqs_av_qual_cutoff";
        } else {
            next; # discard
        }
    }


    if ( $cons_nqs_mm_rate > $cons_av_nqs_mm_cutoff ) {
        # filter indel out: consensus nqs window too messy (probably "strange" indel)
        if ( $mode == 0 ) { # annotate
            $annotation .= "\tAUTOFILTER_CONS_NQS_MM_$cons_av_nqs_mm_cutoff";
        } else {
            next; # discard
        }
    }

    print "$_$annotation\n";


}


close $input_stream if ( $calls ne "-" );

