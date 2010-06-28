#!/usr/bin/perl -w

use strict;
use lib "$ENV{'STING_DIR'}/perl";
use StingArgs;
use DistributedMake;
use DataTable;
use Data::Dumper;

my %args = &getCommandArguments(
    "OUT_DIR" => undef,
    "SEQUENCE_INDEX" => "/broad/1KG/DCC/ftp/sequence.index",
    "ALIGNMENT_INDEX" => "/broad/1KG/DCC/ftp/alignment.index"
);

my @sequenceIndex = &readTable($args{'SEQUENCE_INDEX'}, "header" => 1, "delimiter" => "\t");
my @alignmentIndex = &readTable($args{'ALIGNMENT_INDEX'}, "header" => 1, "delimiter" => "\t");

my %results;

foreach my $sequenceEntry (@sequenceIndex) {
    my %sequenceEntry = %{$sequenceEntry};

    if ($sequenceEntry{'STUDY_NAME'} =~ /1000/ && $sequenceEntry{'ANALYSIS_GROUP'} =~ /low coverage/ && $sequenceEntry{'FASTQ_FILE'} !~ /_\d.filt.fastq.gz/) {
        ${${$results{'all'}}{'sequencer'}}{$sequenceEntry{'INSTRUMENT_PLATFORM'}}++;
        ${${$results{$sequenceEntry{'STUDY_ID'}}}{'sequencer'}}{$sequenceEntry{'INSTRUMENT_PLATFORM'}}++;

        ${$results{'all'}}{'population'} = "all";
        ${$results{$sequenceEntry{'STUDY_ID'}}}{'population'} = $sequenceEntry{'POPULATION'};

        ${$results{'all'}}{'used_lanes'} += ($sequenceEntry{'WITHDRAWN'} == 0);
        ${$results{$sequenceEntry{'STUDY_ID'}}}{'used_lanes'} += ($sequenceEntry{'WITHDRAWN'} == 0);

        ${$results{'all'}}{'unused_lanes'} += ($sequenceEntry{'WITHDRAWN'} != 0);
        ${$results{$sequenceEntry{'STUDY_ID'}}}{'unused_lanes'} += ($sequenceEntry{'WITHDRAWN'} != 0);

        if ($sequenceEntry{'WITHDRAWN'} == 0) {
            ${${${$results{'all'}}{'lanes_per_samples'}}}{$sequenceEntry{'SAMPLE_NAME'}}++;
            ${${${$results{$sequenceEntry{'STUDY_ID'}}}{'lanes_per_samples'}}}{$sequenceEntry{'SAMPLE_NAME'}}++;
        } else {
            ${${$results{'all'}}{'withdrawn_comment'}}{$sequenceEntry{'COMMENT'}}++;
            ${${$results{$sequenceEntry{'STUDY_ID'}}}{'withdrawn_comment'}}{$sequenceEntry{'COMMENT'}}++;
        }

        ${${$results{'all'}}{'lane_parities'}}{$sequenceEntry{'LIBRARY_LAYOUT'}}++;
        ${${$results{$sequenceEntry{'STUDY_ID'}}}{'lane_parities'}}{$sequenceEntry{'LIBRARY_LAYOUT'}}++;

        ${${$results{'all'}}{'submission_date'}}{$sequenceEntry{'SUBMISSION_DATE'}}++;
        ${${$results{$sequenceEntry{'STUDY_ID'}}}{'submission_date'}}{$sequenceEntry{'SUBMISSION_DATE'}}++;

        ${${$results{'all'}}{'center_names'}}{$sequenceEntry{'CENTER_NAME'}}++;
        ${${$results{$sequenceEntry{'STUDY_ID'}}}{'center_names'}}{$sequenceEntry{'CENTER_NAME'}}++;
    }
}

foreach my $studyid (keys(%results)) {
    my %study = %{${${$results{$studyid}}{'lanes_per_samples'}}};
    my @samples = keys(%study);
    my $numSamples = scalar(@samples);

    ${$results{$studyid}}{'num_samples'} = $numSamples;
}

print Dumper(\%results);
