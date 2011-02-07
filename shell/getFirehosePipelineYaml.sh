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

FIREHOSE_ANNOTATIONS=(reference_file interval_list \
  sample_id recalibrated_bam_file squid_project collaborator_id)

# YAML templates

# Project YAML template, once per file.
PROJECT_YAML_TEMPLATE='"\n\
  project: {\n\
    name: '"$ENTITY_SET_ID"',\n\
    referenceFile: %s,\n\
    genotypeDbsnpFile: %s,\n\
    evalDbsnpFile: %s,\n\
    refseqTable: %s,\n\
    intervalList: %s\n\
  },", $1, genotypeDbsnp, evalDbsnp, refseq, $2'

# Project YAML template, once per sample.
SAMPLE_YAML_TEMPLATE='"\n\
    {\n\
      id: %s,\n\
      bamFiles: { cleaned: %s },\n\
      tags: {\n\
        SQUIDProject: %s,\n\
        CollaboratorID: %s\n\
      }\n\
    }", $3, $4, $5, $6'

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
    refseq_dir = "/humgen/gsa-hpprojects/GATK/data/Annotations/refseq/";
    dbsnp_dir = "/humgen/gsa-hpprojects/GATK/data/";

    # add hg18 specific files to awk associative arrays
    genotypeDbsnps["Homo_sapiens_assembly18.fasta"] = dbsnp_dir "dbsnp_129_hg18.rod";
    evalDbsnps["Homo_sapiens_assembly18.fasta"] = dbsnp_dir "dbsnp_129_hg18.rod";
    refseqs["Homo_sapiens_assembly18.fasta"] = refseq_dir "refGene-big-table-hg18.txt";

    # add hg19 specific files to awk associative arrays
    genotypeDbsnps["Homo_sapiens_assembly19.fasta"] = dbsnp_dir "dbsnp_132_b37.vcf";
    evalDbsnps["Homo_sapiens_assembly19.fasta"] = dbsnp_dir "dbsnp_129_b37.rod";
    refseqs["Homo_sapiens_assembly19.fasta"] = refseq_dir "refGene-big-table-hg19.txt";

    printf "{"
}
{
    if (NR == 1) {
        # Based on the reference of the first sample, specify the dbsnps and refseq tables.

        reference_part_count = split($1, reference_parts, "/")
        reference_name = reference_parts[reference_part_count];

        genotypeDbsnp = genotypeDbsnps[reference_name];
        evalDbsnp = evalDbsnps[reference_name];
        refseq = refseqs[reference_name];

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

#hg19=`grep "assembly19" -c $PIPELINE_YAML_FILE`

# NOTE: DBSNP's are populated via AWK's BEGIN block above.
#if [ "$hg19" -ne 0 ]; then
#    sed 's/\/humgen.*rod/\/humgen\/gsa-hpprojects\/GATK\/data\/dbsnp_132_b37.vcf/' $PIPELINE_YAML_FILE > yaml2
#    mv yaml2 $PIPELINE_YAML_FILE
#fi

# NOTE: Renamed "recalibrated" to "cleaned" in SAMPLE_YAML_TEMPLATE above.
#sed 's/recalibrat/clean/' $PIPELINE_YAML_FILE > yaml2
#mv yaml2 $PIPELINE_YAML_FILE
