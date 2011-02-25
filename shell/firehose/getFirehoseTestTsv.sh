#!/bin/sh

# Downloads a set of samples from Firehose using the obsolete Firehose Test Harness and generate a TSV file with the outputs.

ENTITY_SET_ID=$1
ENTITY_SET_TYPE=Sample_Set
ENTITY_TYPE=Sample

if [ "$ENTITY_SET_ID" == "" ]; then
    echo "Usage: $0 <Sample_Set_Name>" >&2
    exit 1
fi

# Firehose variables

FIREHOSE_SOURCE_HOME=/humgen/gsa-firehose/firehose/source
CGA_HOME=$FIREHOSE_SOURCE_HOME/CancerGenomeAnalysis
FIREHOSE_TEST_HARNESS="python $CGA_HOME/analysis_pipeline/scripts/firehose_test_harness.py"
FIREHOSE_HOST=firehose
FIREHOSE_PORT=8080
FIREHOSE_DOMAIN=gsa
FIREHOSE_WORKSPACE=trunk

# TSV file to write

PIPELINE_TSV_FILE=$ENTITY_SET_ID.tsv

# Annotations to pull down from Firehose

FIREHOSE_ANNOTATIONS=(reference_file interval_list sample_id recalibrated_bam_file squid_project collaborator_id)

index=0
count=${#FIREHOSE_ANNOTATIONS[@]}
TSV_HEADER=""
FIREHOSE_VARIABLES=""
TAB='	'

# Build the tab separated list of firehose arguments

while [ "$index" -lt "$count" ]; do
    if [ "$FIREHOSE_VARIABLES" != "" ]; then
        FIREHOSE_VARIABLES=$FIREHOSE_VARIABLES$TAB
        TSV_HEADER=$TSV_HEADER$TAB
    fi
    FIREHOSE_VARIABLES=$FIREHOSE_VARIABLES'${'${FIREHOSE_ANNOTATIONS[$index]}'}'
    TSV_HEADER=$TSV_HEADER${FIREHOSE_ANNOTATIONS[$index]}
    let "index = $index + 1"
done

# Retrieve all the required variables via the test harness.
$FIREHOSE_TEST_HARNESS \
    -d $FIREHOSE_DOMAIN -w $FIREHOSE_WORKSPACE \
    -t $ENTITY_TYPE -f $ENTITY_SET_ID -y $ENTITY_SET_TYPE \
    "echo '$FIREHOSE_VARIABLES'" && \
\
# Generate tsv header
echo "$TSV_HEADER" > $PIPELINE_TSV_FILE \
# Generate tsv from firehose output
. firehose-populated-commands.sh >> $PIPELINE_TSV_FILE

EXIT_CODE=$?

if [[ $EXIT_CODE -ne 0 ]]; then
    echo "" >&2
    echo "The Firehose test harness failed with exit code:" $EXIT_CODE >&2
    echo 'Check the name of your Sample_Set or try using the newer getFirehoseCurlTsv.sh' >&2
    exit $EXIT_CODE
fi
