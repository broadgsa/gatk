#!/usr/bin/env python

# Executable files.  Please update these to match the installed locations of your tools.
samtools_exe='/seq/dirseq/samtools/current/samtools'
java_exe='/broad/tools/Linux/x86_64/pkgs/jdk_1.6.0_12/bin/java'
R_exe='/util/bin/R'
GATK_jar='/humgen/gsa-scr1/hanna/src/StingWorking/dist/GenomeAnalysisTK.jar'

# Location of the resource files distributed with the recalibration tool.
resources='resources/'

# Where does the reference live?
reference_base = resources + 'Homo_sapiens_assembly18'
reference      = reference_base + '.fasta'
reference_dict = reference_base + '.dict'
reference_fai  = reference_base + '.fasta.fai'

# Where does DBSNP live?
dbsnp = resources + 'dbsnp.rod.out'

import sys,os

def exit(msg,errorcode):
    print msg
    sys.exit(errorcode)

def check_input_file_available(filename,description):
    if not os.access(filename,os.R_OK):
        exit('Unable to access %s %s' % (description,filename),1)

if len(sys.argv) < 3:
    exit('Usage: python RecalQual.py <input bam file> <calibrated output bam file>',1)

# check that the input bam file exists, and that the bam is indexed.
bam = sys.argv[1]
bam_index = bam + '.bai'

check_input_file_available(bam,'reads file')
check_input_file_available(bam_index,'reads index file')

# parse the user's calibration output file requirements
calibrated_bam = sys.argv[2]
calibrated_bam_index = calibrated_bam + '.bai'

# check that the fasta and supporting files are available
check_input_file_available(reference,'reference file')
check_input_file_available(reference_dict,'reference dictionary')
check_input_file_available(reference_fai,'reference index file')

# check that the dbsnp is available
check_input_file_available(dbsnp,'dbsnp file')

# make an output directory for temporary files
if not os.path.isdir('output'):
    os.mkdir('output')

# assemble the required program arguments
gatk_base_cmdline = ' '.join((java_exe,'-ea','-jar',GATK_jar,'-R',reference,'--DBSNP',dbsnp,'-l INFO','-L chrM'))

generate_covariates         = ' '.join((gatk_base_cmdline,'-T CountCovariates','-I',bam,'-mqs 40','--OUTPUT_FILEROOT output/output','--CREATE_TRAINING_DATA','--MIN_MAPPING_QUALITY 1'))
compute_logistic_regression = ' '.join(('python','LogisticRegressionByReadGroup.py','output/output.covariate_counts.csv','output/linear_regression_results.out'))
apply_logistic_regression   = ' '.join((gatk_base_cmdline,'-T LogisticRecalibration','-I',bam,'-logisticParams output/linear_regression_results.out','-outputBAM',calibrated_bam))
index_calibrated_bamfile    = ' '.join((samtools_exe,'index',calibrated_bam))

os.system(generate_covariates)
os.system(compute_logistic_regression)
os.system(apply_logistic_regression)
os.system(index_calibrated_bamfile)
