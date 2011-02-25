#!/bin/sh

# Downloads a set of samples from Firehose and generates a YAML file.

DIR=`dirname $0`
if [ "$2" == "" ]; then
    $DIR/getFirehoseTestTsv.sh $1 && $DIR/pipelineTsvToYaml.sh $1.tsv
else
    $DIR/getFirehoseCurlTsv.sh $1 $2 && $DIR/pipelineTsvToYaml.sh $1.tsv
fi
