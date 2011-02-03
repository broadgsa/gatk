#!/bin/sh

SAMPLE_SET=$1
STING_HOME=/humgen/gsa-pipeline/.repository
REFSEQ_TABLE=/humgen/gsa-hpprojects/GATK/data/Annotations/refseq/refGene-big-table-hg19.txt
JOB_QUEUE=week
SHORT_QUEUE=hour
TMP_DIR=$PWD/tmp

mkdir $PWD/tmp
source /broad/software/scripts/useuse
use LSF
use R-2.10
use Oracle-full-client
use .cx-oracle-5.0.2-python-2.6.5-oracle-full-client-11.1

java -Djava.io.tmpdir="$TMP_DIR" -jar "$STING_HOME"/dist/Queue.jar -jobQueue "$JOB_QUEUE" -shortJobQueue "$SHORT_QUEUE" -jobProject "$SAMPLE_SET" -jobPrefix "$SAMPLE_SET" -tearScript "$STING_HOME"/R/DataProcessingReport/GetTearsheetStats.R -S "$STING_HOME"/scala/qscript/playground/fullCallingPipeline.q -Y "$SAMPLE_SET".yaml -refseqTable "$REFSEQ_TABLE" --gatkjar "$STING_HOME"/dist/GenomeAnalysisTK.jar -log queue_log.txt -statusTo corin -bsub $2
