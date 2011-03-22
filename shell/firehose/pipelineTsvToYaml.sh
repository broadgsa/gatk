#!/bin/sh

# Uses awk to generate a YAML file from a TSV.

# In the awk script and templates:
#   - Variables starting with a '$' are columns in the TSV
#   - Variables without a '$' are pre-calculated from the first row of data


# TSV file to read
PIPELINE_TSV_FILE=$1

if [ "$PIPELINE_TSV_FILE" == "" ]; then
    echo "Usage: $0 <Set_Name>.tsv" >&2
    exit 1
fi

ROW_COUNT=(`wc -l $PIPELINE_TSV_FILE`)
if [[ ${ROW_COUNT[0]} -lt 2 ]]; then
    echo "Header plus data not found in tsv: $PIPELINE_TSV_FILE" >&2
    exit 1
fi

# YAML file to write
PIPELINE_YAML_FILE=${PIPELINE_TSV_FILE%.tsv}.yaml

# YAML templates

# Project YAML template, once per file.
PROJECT_YAML_TEMPLATE='"\n\
  project: {\n\
    name: %s,\n\
    referenceFile: %s,\n\
    genotypeDbsnp: %s,\n\
    evalDbsnp: %s,\n\
    refseqTable: %s,\n\
    intervalList: %s\n\
  },", projectName, $referenceFile, genotypeDbsnp, evalDbsnp, refseq, $intervalList'

# Project YAML template, once per sample.
SAMPLE_YAML_TEMPLATE='"\n\
    {\n\
      id: %s,\n\
      bamFiles: { cleaned: %s },\n\
      tags: {\n\
        SQUIDProject: %s,\n\
        CollaboratorID: %s\n\
      }\n\
    }", $sampleId, $bamFile, $squidProject, $collaboratorId'

TEST_AWK_COUNT=`echo '\n' | awk '{print $0}' | wc -c`
if [ "$TEST_AWK_COUNT" -eq 2 ]; then
    # Strip the extra \n from the lines if awk of \n is
    # a newline and not the two characters slash-n (on mac)
    PROJECT_YAML_TEMPLATE="${PROJECT_YAML_TEMPLATE//\\\n/}"
    SAMPLE_YAML_TEMPLATE="${SAMPLE_YAML_TEMPLATE//\\\n/}"
fi

# Generate yaml from tsv
awk '
{
    if (NR == 1) {
        tsvFile = "'$PIPELINE_TSV_FILE'"

        # Set the project name to the TSV file minus the directory and the .tsv
        projectName = tsvFile
        sub(/\/.*\//, "", projectName)
        sub(/\.tsv/, "", projectName)

        # Read the column headers and figure out the index of each column name.
        for (i=1; i<=NF; i++)
            columnFields[tolower($i)] = i

        referenceFile = columnFields["reference_file"]
        intervalList = columnFields["interval_list"]
        sampleId = columnFields["sample_id"]
        squidProject = columnFields["squid_project"]
        collaboratorId = columnFields["collaborator_id"]

        for (key in columnFields)
            if (key ~ "bam_file")
                bamFile = columnFields[key]

        if (referenceFile == "") {
            print "ERROR: Column header reference_file missing from " tsvFile > "/dev/stderr"
            exitWithError = 1
        }

        if (intervalList == "") {
            print "ERROR: Column header interval_list missing from " tsvFile > "/dev/stderr"
            exitWithError = 1
        }

        if (sampleId == "") {
            print "ERROR: Column header sample_id missing from " tsvFile > "/dev/stderr"
            exitWithError = 1
        }

        if (squidProject == "") {
            print "ERROR: Column header squid_project missing from " tsvFile > "/dev/stderr"
            exitWithError = 1
        }

        if (collaboratorId == "") {
            print "ERROR: Column header collaborator_id missing from " tsvFile > "/dev/stderr"
            exitWithError = 1
        }

        if (bamFile == "") {
            print "ERROR: Column header *bam_file* missing from " tsvFile > "/dev/stderr"
            exitWithError = 1
        }

        if (exitWithError) {
            exit 1
        }

        refseqDir = "/humgen/gsa-hpprojects/GATK/data/Annotations/refseq/"
        dbsnpDir = "/humgen/gsa-hpprojects/GATK/data/"

        # add hg18 specific files to awk associative arrays
        genotypeDbsnps["Homo_sapiens_assembly18.fasta"] = dbsnpDir "dbsnp_129_hg18.rod"
        evalDbsnps["Homo_sapiens_assembly18.fasta"] = dbsnpDir "dbsnp_129_hg18.rod"
        refseqs["Homo_sapiens_assembly18.fasta"] = refseqDir "refGene-big-table-hg18.txt"

        # add hg19 specific files to awk associative arrays
        genotypeDbsnps["Homo_sapiens_assembly19.fasta"] = dbsnpDir "dbsnp_132_b37.vcf"
        evalDbsnps["Homo_sapiens_assembly19.fasta"] = dbsnpDir "dbsnp_129_b37.vcf"
        refseqs["Homo_sapiens_assembly19.fasta"] = refseqDir "refGene-big-table-hg19.txt"

        printf "{"
    } else {
        missingValue = 0
        if ($referenceFile == "") missingValue = 1
        if ($intervalList == "") missingValue = 1
        if ($sampleId == "") missingValue = 1
        if ($squidProject == "") missingValue = 1
        if ($collaboratorId == "") missingValue = 1
        if ($bamFile == "") missingValue = 1

        if (missingValue) {
            print "WARNING: Skipping row which does not have all values: " $0 > "/dev/stderr"
        } else {
            if (NR == 2) {
                # Based on the reference of the first sample, specify the dbsnps and refseq tables.

                referencePartCount = split($referenceFile, referenceParts, "/")
                referenceName = referenceParts[referencePartCount]

                genotypeDbsnp = genotypeDbsnps[referenceName]
                evalDbsnp = evalDbsnps[referenceName]
                refseq = refseqs[referenceName]

                printf '"$PROJECT_YAML_TEMPLATE"'
                printf "\n  samples: ["
            } else {
                printf ","
            }
            printf '"$SAMPLE_YAML_TEMPLATE"'
        }
    }
}
END {
    if (NR > 0)
        printf "\n  ]"
    print "\n}"
}' "$PIPELINE_TSV_FILE" > "$PIPELINE_YAML_FILE"
