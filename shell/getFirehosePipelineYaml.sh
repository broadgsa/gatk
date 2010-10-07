#!/bin/sh

# Downloads a set of samples from Firehose using the Firehose Test Harness and awk to generate a YAML file.

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

# YAML file to write

PIPELINE_YAML_FILE=$ENTITY_SET_ID.yaml

# Annotations to pull down from Firehose

FIREHOSE_ANNOTATIONS=(reference_file dbsnp_file interval_list \
  sample_id recalibrated_bam_file squid_project collaborator_id)

# YAML templates

PROJECT_YAML_TEMPLATE='"\n\
  project: {\n\
    name: '"$ENTITY_SET_ID"',\n\
    referenceFile: %s,\n\
    dbsnpFile: %s,\n\
    intervalList: %s\n\
  },", $1, $2, $3'

SAMPLE_YAML_TEMPLATE='"\n\
    {\n\
      id: %s,\n\
      bamFiles: { recalibrated: %s },\n\
      tags: {\n\
        SQUIDProject: %s,\n\
        CollaboratorID: %s\n\
      }\n\
    }", $4, $5, $6, $7'

TEST_AWK_COUNT=`echo '\n' | awk '{print $0}' | wc -c`
if [ "$TEST_AWK_COUNT" -eq 2 ]; then
    # Strip the extra \n from the lines if awk of \n is
    # a newline and not the two characters slash-n (on mac)
    PROJECT_YAML_TEMPLATE="${PROJECT_YAML_TEMPLATE//\\\n/}"
    SAMPLE_YAML_TEMPLATE="${SAMPLE_YAML_TEMPLATE//\\\n/}"
fi

index=0
count=${#FIREHOSE_ANNOTATIONS[@]}
FIREHOSE_VARIABLES=""
TAB='	'

# Build the tab separated list of firehose arguments

while [ "$index" -lt "$count" ]; do
    if [ "$FIREHOSE_VARIABLES" != "" ]; then
        FIREHOSE_VARIABLES=$FIREHOSE_VARIABLES$TAB
    fi
    FIREHOSE_VARIABLES=$FIREHOSE_VARIABLES'${'${FIREHOSE_ANNOTATIONS[$index]}'}'
    let "index = $index + 1"
done

# Retrieve all the required variables and run the pipeline in Queue.
$FIREHOSE_TEST_HARNESS \
    -d $FIREHOSE_DOMAIN -w $FIREHOSE_WORKSPACE \
    -t $ENTITY_TYPE -f $ENTITY_SET_ID -y $ENTITY_SET_TYPE \
    "echo '$FIREHOSE_VARIABLES'" && \
\
# Generate yaml from firehose output
. firehose-populated-commands.sh | awk '
BEGIN {
    printf "{"
}
{
    if (NR == 1) {
      printf '"$PROJECT_YAML_TEMPLATE"'
      printf "\n  samples: ["
    } else {
        printf ","
    }
    printf '"$SAMPLE_YAML_TEMPLATE"'
}
END {
    if (NR > 0)
        printf "\n  ]"
    print "\n}"
}' > $PIPELINE_YAML_FILE
