#!/bin/sh

SAMPLE_SET=$1
STING_HOME=/humgen/gsa-firehose2/pipeline/repositories/StingProduction
REFSEQ_TABLE=/humgen/gsa-hpprojects/GATK/data/Annotations/refseq/refGene-big-table-hg19.txt
JOB_QUEUE=gsa
SHORT_QUEUE=gsa
BIG_QUEUE=gsa
TMP_DIR=$PWD/tmp

mkdir $PWD/tmp
source /broad/software/scripts/useuse
use LSF
use R-2.10
use Oracle-full-client
use .cx-oracle-5.0.2-python-2.6.5-oracle-full-client-11.1

java -Djava.io.tmpdir="$TMP_DIR" -jar "$STING_HOME"/dist/Queue.jar -jobQueue "$JOB_QUEUE" -shortJobQueue "$SHORT_QUEUE" -bigMemQueue "$BIG_QUEUE" -jobProject "$SAMPLE_SET" -jobPrefix "$SAMPLE_SET" -tearScript ~/Sting/R/DataProcessingReport/GetTearsheetStats.R -S ~/Sting/scala/qscript/fullCallingPipeline.q -Y "$SAMPLE_SET".yaml -refseqTable "$REFSEQ_TABLE" --gatkjar "$STING_HOME"/dist/GenomeAnalysisTK.jar -log queue_log.txt -statusTo corin -bsub $2 
