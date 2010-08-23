#!/bin/sh

ENTITY_TYPE=$1
ENTITY_ID=$2
FIREHOSE_TOKEN=$3

FIREHOSE_WORKSPACE=trunk
FIREHOSE_HOST=firehose
FIREHOSE_PORT=8080
FIREHOSE_DOMAIN=gsa
CLEAN_BAM_ANNOTATION=clean_bam_file
DEFAULT_QUEUE=gsa
SHORT_QUEUE=hour
FIREHOSE_SOURCE_HOME=/humgen/gsa-firehose/firehose/source
STING_HOME=$FIREHOSE_SOURCE_HOME/Sting
CGA_HOME=$FIREHOSE_SOURCE_HOME/CancerGenomeAnalysis
RUN=-run
TMP_DIR=/broad/shptmp/$USER
JOB_QUEUE=gsa
INTERVAL_COUNT=15
REALIGNER_COUNT=50
QUEUE_JAR=$STING_HOME/dist/Queue.jar
GATK_JAR=$STING_HOME/dist/GenomeAnalysisTK.jar
FIX_MATES_JAR=/seq/software/picard/1.194/bin/FixMateInformation.jar
SAMTOOLS=/seq/dirseq/samtools/samtools-0.1.7-5/samtools
FIREHOSE_IMPORT_JAR=$CGA_HOME/analysis_pipeline/tools/dist/ImportSingleValue.jar

CLEAN_BAM_FILE_SCRIPT=$STING_HOME/scala/qscript/kshakir/CleanBamFile.scala
FIREHOSE_TEST_HARNESS="python $CGA_HOME/analysis_pipeline/scripts/firehose_test_harness.py"
SPLIT_BY_CONTIG="python $STING_HOME/python/splitIntervalsByContig.py"
SPLIT_BY_INTERVALS=$STING_HOME/shell/splitIntervals.sh
MERGE_TEXT=$STING_HOME/shell/mergeText.sh
CONVERT_BLACKLIST=$CGA_HOME/analysis_pipeline/genepattern/modules/ConvertDependentsList/convert_dependents_list.sh

# Record the Sting version number.
svnversion $STING_HOME > stingversion.txt
# If Sting has been modified, record the differences.
grep -E '^[0-9]+$' stingversion.txt || svn diff $STING_HOME > stingdiff.txt

# Record the CGA version number.
svnversion $CGA_HOME > cgaversion.txt
# If CGA has been modified, record the differences.
grep -E '^[0-9]+$' cgaversion.txt || svn diff $CGA_HOME > cgadiff.txt

# Try to retrieve the blacklist.  If it fails then set it to "".
$FIREHOSE_TEST_HARNESS -d $FIREHOSE_DOMAIN -w $FIREHOSE_WORKSPACE -n $ENTITY_ID -t $ENTITY_TYPE 'BLACKLIST="${read_group_blacklist}"' && . firehose-populated-commands.sh || BLACKLIST=""

# Retrieve all the required variables and run the pipeline in Queue.
$FIREHOSE_TEST_HARNESS -d $FIREHOSE_DOMAIN -w $FIREHOSE_WORKSPACE -n $ENTITY_ID -t $ENTITY_TYPE 'REFERENCE_FILE="${reference_file}";BAM_FILE="${recalibrated_bam_file}";DBSNP_FILE="${dbsnp_file}";INTERVAL_FILE="${interval_list}";DATABASE_ID="${database_id}"' && . firehose-populated-commands.sh && \
\
JOB_PREFIX=Q-$ENTITY_ID && \
\
java \
-Djava.io.tmpdir="$TMP_DIR" \
-jar "$QUEUE_JAR" \
-S "$CLEAN_BAM_FILE_SCRIPT" \
-bsub -bsubWait -skipUpToDate \
-jobQueue "$JOB_QUEUE" \
-jobPrefix "$JOB_PREFIX" \
-gatk "$GATK_JAR" \
-base "$ENTITY_ID" \
-R "$REFERENCE_FILE" \
-I "$BAM_FILE" \
-D "$DBSNP_FILE" \
-L "$INTERVAL_FILE" \
-RTCSS "$SPLIT_BY_CONTIG" \
-RTCSC "$INTERVAL_COUNT" \
-IRSS "$SPLIT_BY_INTERVALS" \
-IRSC "$REALIGNER_COUNT" \
-MTS "$MERGE_TEXT" \
-fixMates "$FIX_MATES_JAR" \
-samtools "$SAMTOOLS" \
-RGBLS "$CONVERT_BLACKLIST" \
-RGBL "$BLACKLIST" \
-importJar "$FIREHOSE_IMPORT_JAR" \
-shortQueue "$SHORT_QUEUE" \
-FHHost "$FIREHOSE_HOST" \
-FHPort "$FIREHOSE_PORT" \
-FHDomain "$FIREHOSE_DOMAIN" \
-FHToken "$FIREHOSE_TOKEN" \
-bamFHEType "$ENTITY_TYPE" \
-bamFHEID "$DATABASE_ID" \
-bamFHAnn "$CLEAN_BAM_ANNOTATION" \
$RUN || ( echo Job failed. Check stdout.txt for more info. >&2; exit 1; )
