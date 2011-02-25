#!/bin/sh

# Downloads a set of samples from Firehose using the Firehose API and generate a TSV file with the outputs.
# see: http://iwww.broadinstitute.org/cancer/cga/wiki/index.php/GetAnnotations

ENTITY_SET_ID=$1
ENTITY_SET_TYPE=Sample_Set
ENTITY_TYPE=Sample
PASSWORD_FILE=$2

if [ "$ENTITY_SET_ID" == "" ]; then
    EXIT_USAGE=1
fi

if [ "$PASSWORD_FILE" == "" ]; then
    echo 'Missing password file with the contents: "-u <user>:<pass>"' >&2
    EXIT_USAGE=1
fi

if [ $EXIT_USAGE ]; then
    echo "Usage: $0 <Sample_Set_Name> <Curl_Password_File>" >&2
    exit 1
fi

# Firehose variables

FIREHOSE_HOST=firehose
FIREHOSE_PORT=8080
FIREHOSE_DOMAIN=gsa
FIREHOSE_WORKSPACE=trunk

# TSV file to write

PIPELINE_TSV_FILE=$ENTITY_SET_ID.tsv

# Annotations to pull down from Firehose

FIREHOSE_ANNOTATIONS=(reference_file interval_list recalibrated_bam_file squid_project collaborator_id)

index=0
count=${#FIREHOSE_ANNOTATIONS[@]}
FIREHOSE_VARIABLES=""

# Build the tab separated list of firehose arguments

while [ "$index" -lt "$count" ]; do
    FIREHOSE_VARIABLES=$FIREHOSE_VARIABLES'&annotationTypes='${FIREHOSE_ANNOTATIONS[$index]}
    let "index = $index + 1"
done

curl --fail -sL -K "$PASSWORD_FILE" -o "$PIPELINE_TSV_FILE" \
    "http://$FIREHOSE_HOST:$FIREHOSE_PORT/$FIREHOSE_DOMAIN/ws/entity/getAnnotations/$ENTITY_TYPE?entityNames=$ENTITY_SET_ID&filterSetType=$ENTITY_SET_TYPE&workspaceName=$FIREHOSE_WORKSPACE$FIREHOSE_VARIABLES" || \

EXIT_CODE=$?

if [[ $EXIT_CODE -ne 0 ]]; then
    echo "curl failed with exit code:" $EXIT_CODE >&2
    echo 'Check the name of your Sample_Set and that your password file '$PASSWORD_FILE' is setup correctly with: "-u <user>:<pass>"' >&2
    echo "If that doesn't work make sure you can login to the firehose website: http://$FIREHOSE_HOST:$FIREHOSE_PORT/$FIREHOSE_DOMAIN" >&2
    exit $EXIT_CODE
fi
