#!/usr/bin/env python

# Executable files.  Please update these to match the installed locations of your tools.
samtools_exe='/seq/dirseq/samtools/current/samtools'
java_exe='/broad/tools/Linux/x86_64/pkgs/jdk_1.6.0_12/bin/java'
R_exe="/broad/tools/apps/R-2.6.0/bin/Rscript"

# Any special site-specific arguments to pass the JVM.
jvm_args='-ea'

# Which platforms should the calibration tool be run over?
platforms=['illumina']

# Where to put the output created as part of recalibration.  
# If editing, please end this variable with a trailing slash.
output_root = './'

# Location of the resource files distributed with the recalibration tool.  
# If editing, please end this variable with a trailing slash.
resources='recalibratorResources/'
#resources='/broad/1KG/DCC_recal/ReadQualityRecalibrator/resources/'

import getopt,glob,os,sys
import LogisticRegressionByReadGroup

# if the provided resources path is relative, convert it to a path
# relative to the calling script rather than the cwd.
if not os.path.isabs(resources):
    application_path = os.path.dirname(sys.argv[0])
    resources = os.path.abspath(application_path) + '/' + resources

# Where does the reference live?
reference_base = resources + 'Homo_sapiens_assembly18'
reference      = reference_base + '.fasta'
reference_dict = reference_base + '.dict'
reference_fai  = reference_base + '.fasta.fai'

# Where does DBSNP live?
dbsnp = resources + 'dbsnp_129_hg18.rod'

# Where are the application files required to run the recalibration?
#gatk = resources + 'gatk/GenomeAnalysisTK.jar'
gatk = '/home/radon01/depristo/dev/GenomeAnalysisTK_stable/trunk/dist/GenomeAnalysisTK.jar'

logistic_regression_script = resources + 'logistic_regression.R'
empirical_vs_reported_grapher = resources + 'plot_q_emp_stated_hst.R'

# Assemble the platform list into command-line arguments.
platform_args = ' '.join(['-pl %s' % platform for platform in platforms])

def output_file_needs_update(output_file, input_files = []):
    if reprocessFiles:                          # we are always reprocessing files 
        return True
    if not os.path.exists(output_file):         # the file doesn't exist
        return True
    else: 
        return False 

def executeCommand(cmd):
    if not dryRun:
        return os.system(cmd)
    else:
        print '  => Would execute', cmd
        return 0

def exit(msg,errorcode):
    print msg
    sys.exit(errorcode)

def check_input_file_available(filename,description):
    if not os.access(filename,os.R_OK):
        exit('Unable to access %s %s' % (description,filename),1)

def graph_file(graph_script,graph_data):
    'Graph the given data using the given script.  Leave the data in the output directory.'
    check_input_file_available(graph_script,'%s R graphing script' % graph_script)
    check_input_file_available(graph_data,'%s graphing data' % graph_data)
    result = executeCommand(' '.join((R_exe,graph_script,graph_data)))
    if result != 0:
        exit('Unable to graph data: %s' % graph_data,1)

def recalibrate():
    'Recalibrate the given bam file'
    # generate the covariates
    print 'generating covariates'
    linear_regression_results = output_dir + 'linear_regression_results.out'
    if output_file_needs_update(linear_regression_results):
        generate_covariates = ' '.join((gatk_base_cmdline,'-T CountCovariates','-I',bam,'-mqs 40','--OUTPUT_FILEROOT',output_dir+'initial','--CREATE_TRAINING_DATA','--MIN_MAPPING_QUALITY 1',platform_args))
        returncode = executeCommand(generate_covariates)
        if returncode != 0:
            exit('Unable to generate covariates',1)
    
        # compute the logistic regression
        if not dryRun:
            print 'computing the logistic regression'
            LogisticRegressionByReadGroup.compute_logistic_regression(output_dir + 'initial.covariate_counts.csv',linear_regression_results,R_exe,logistic_regression_script)
    else:
        print 'Logistic recalibration data files already generated', linear_regression_results
  
    if output_file_needs_update(calibrated_bam):
        # apply the logistic regression, writing the output data to calibrated_bam
        print 'applying the correction to the reads'
        apply_logistic_regression = ' '.join((gatk_base_cmdline,'-T LogisticRecalibration','-I',bam,'-logisticParams',output_dir+'linear_regression_results.out','-outputBAM',calibrated_bam))
        returncode = executeCommand(apply_logistic_regression)
        if returncode != 0:
            exit('Unable to apply logistic regression',1)
    else:
        print 'Recalibrated BAM file already generated', calibrated_bam

    # index the calibrated bam
    if output_file_needs_update(calibrated_bam + '.bai'):
        print 'indexing the calibrated bam'
        index_calibrated_bamfile = ' '.join((samtools_exe,'index',calibrated_bam))
        returncode = executeCommand(index_calibrated_bamfile)
        if returncode != 0:
            exit('Unable to index calibrated bamfile',1)
    else:
        print 'Recalibrated BAM index file already generated', calibrated_bam + '.bai'

    print 'Recalibration complete!  Calibrated bam is available here: ' + calibrated_bam    

def evaluate():
    'Evaluate recalibration results.'
    print 'Evaluating recalibration results'

    recalibrated_regression_results = output_dir + 'recalibrated' + 'linear_regression_results.out'
    if output_file_needs_update(recalibrated_regression_results):
       # regenerate the covariates
       regenerate_covariates = ' '.join((gatk_base_cmdline,'-T CountCovariates','-I',calibrated_bam,'-mqs 40','--OUTPUT_FILEROOT',output_dir+'recalibrated','--CREATE_TRAINING_DATA','--MIN_MAPPING_QUALITY 1',platform_args))
       print 'regenerating covariates'
       returncode = executeCommand(regenerate_covariates)
       if returncode != 0:
           exit('Unable to regenerate covariates',1)
   
       print 'graphing initial results'
       for filename in glob.glob(output_dir+'initial.*.empirical_v_reported_quality.csv'):
           graph_file(empirical_vs_reported_grapher,filename)
       print 'graphing final results'
       for filename in glob.glob(output_dir+'recalibrated.*.empirical_v_reported_quality.csv'):
           graph_file(empirical_vs_reported_grapher,filename)
    else:
        print 'Evaluated files already generated', recalibrated_regression_results

def usage():
    exit('Usage: python RecalQual.py [--recalibrate] [--evaluate] <input bam file> <calibrated output bam file>',1)

# Try to parse the given command-line arguments.
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],'',['recalibrate','evaluate', 'reprocess', 'dry', 'outputRoot='])
except getopt.GetoptError, err:
    usage()

# An input and an output file is required.  Fail if left unspecified.
if len(args) < 2:
    usage()

# Determine whether to evaluate / recalibrate.
recalibrate_requested = False
evaluate_requested = False
reprocessFiles = False
dryRun = False
for opt,arg in opts:
    if opt == '--recalibrate':
        recalibrate_requested = True
    if opt == '--evaluate':
        evaluate_requested = True
    if opt == '--reprocess':
        reprocessFiles = True
    if opt == '--dry':
        dryRun = True
    if opt == '--outputRoot':
        output_root = arg
        
# Default to 'recalibrate' unless the user specified only evaluate.
do_recalibration = not (evaluate_requested and not recalibrate_requested)
# Only evaluate if the user specifically requested evaluation.
do_evaluation = evaluate_requested

print 'Looking for supplemental files in ' + resources

# check that the input bam file exists, and that the bam is indexed.
bam = args[0]
bam_index = bam + '.bai'

check_input_file_available(bam,'reads file')
check_input_file_available(bam_index,'reads index file')

# parse the user's calibration output file requirements
calibrated_bam = args[1]
calibrated_bam_index = calibrated_bam + '.bai'

# check that the fasta and supporting files are available
check_input_file_available(reference,'reference file')
check_input_file_available(reference_dict,'reference dictionary')
check_input_file_available(reference_fai,'reference index file')

# check that the dbsnp is available
check_input_file_available(dbsnp,'dbsnp file')

# sanity check that the software is available
check_input_file_available(samtools_exe,'samtools')
check_input_file_available(java_exe,'java')
check_input_file_available(R_exe,'R')
check_input_file_available(gatk,'Genome Analysis Toolkit')
check_input_file_available(logistic_regression_script,'logistic regression script')

# make an output directory for temporary files
output_dir=output_root+'output.' + bam[bam.rfind('/')+1:bam.rfind('.bam')] + '/'
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
if not os.path.isdir(output_dir):
    exit('Unable to create output directory ' + output_dir,1)

# assemble the required program arguments
gatk_base_cmdline = ' '.join((java_exe,jvm_args,'-jar',gatk,'-R',reference,'--DBSNP',dbsnp,'-l INFO'))

if do_recalibration:
    recalibrate()

if do_evaluation:
    evaluate()
