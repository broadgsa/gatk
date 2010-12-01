from farm_commands2 import *
import os.path
import sys
from optparse import OptionParser
from datetime import date
import glob
import operator
import faiReader
import math
import shutil
import string
from madPipelineUtils import *

EXCLUDE_CHRS = ['chrM', 'chrY']
EXTRA_GATK_ARGS = ' -XL chrM -XL chrY ' # -XL chrX -XL chrY '
#EXTRA_GATK_ARGS = ' -XL chrM -XL chrY ' # -XL chrX -XL chrY '
VALIDATION_DIR = '/humgen/gsa-hpprojects/GATK/data/Comparisons'
BAM_ROOT = '/humgen/1kg/analysis/bamsForDataProcessingPapers/'
WE_LIST = '/seq/references/HybSelOligos/whole_exome_agilent_designed_120/whole_exome_agilent_designed_120.targets.interval_list'
WGS_FILTER = [['ABFilter', 'AB > 0.75 && DP > 40'], ['DPFilter', 'DP > 120 || SB > -0.10']] # , ['QDFilter', 'QD < 5.0 && DP > 40']]
WE_FILTER = [['ESPStandard', 'AB > 0.75 || QD < 5.0 || HRun > 3 || SB > -0.10']]

#UG_ARGS = "-mbq 20 -mmq 20 -stand_call_conf 50 -stand_emit_conf 10 -hets 0.78e-3 -dcov 10000 -pnrm GRID_SEARCH"
UG_ARGS = "-stand_call_conf 10 -stand_emit_conf 10 --downsampling_type BY_SAMPLE -dcov 1000 -hets 0.78e-3" # experimental arguments for GdA test

class CallTarget:
    def __init__(self, name, bam, interval = '', callArgs = "", b36 = False, optimize = True, filters = [], targetTiTv = 2.07, maxClusters = 16, minQual = 300, tranchToTake = 1):
        self.name = name
        self.bam = bam
        self.interval = interval
        self.callArgs = callArgs
        self.vcfs = []  # list of processed vcf
        self.b36 = b36
        self.filters = filters
        self.optimize = optimize
        self.targetTiTv = targetTiTv
        self.maxClusters = maxClusters
        self.tranchToTake = tranchToTake
        self.minQual = minQual
    
    def getCallArgs(self):
        return self.callArgs # + self.getIntervalArg()
    
    def getIntervalArg(self):
        if self.hasInterval():
            return ' -L ' + self.interval
        else:
            return ''
            
    def hasInterval(self):
        return self.interval != '' and self.interval != None

    def getVcf(self):
        return os.path.join(OPTIONS.dir, self.name + ".vcf")

    def getVcfs(self):
        return self.vcfs

    def addVcf(self, vcf):
        self.vcfs.append(vcf)

    def getBam(self):
        return self.bam

KG_PATH = '/humgen/gsa-hpprojects/1kg/1kg_pilot2/currentBestProjectCalls'
TECH_COMP = '/humgen/gsa-hphome1/kiran/one_off_projects/multiTechComparisons/results/v7/NA12878'

#WGS_INTERVAL = 'chr1'
#WGS_INTERVAL = '-L chr1:1-50,000,000'

def weTarget(name, bam, ignore = '', args = '', filters = None):
    # 3.0 was old target, new is 2.8
    return CallTarget(name, bam, interval = WE_LIST, callArgs = args, filters = WE_FILTER, targetTiTv = 2.8, maxClusters = 12, minQual = 2800, tranchToTake = 10) 

#TARGETS_BY_STRATEGY = [['', ''], ['.OQ', '-OQ'], ['.OQ.noCM', '-OQ -bm THREE_STATE'], ['.noCM', '-bm THREE_STATE']]
TARGETS_BY_STRATEGY = [['', ''], ['.OQ', '-OQ']]
#TARGETS_BY_STRATEGY = [['', '']]
def targetsByStrategy(func, rootName, bam, interval = '', args = '', filters = []):
    def makeTarget(ext, moreArgs):
        name = rootName + ext
        return func(name, bam, interval, args + ' ' + moreArgs, filters = filters)

    if "cleaned" in rootName or "CG" in rootName:
        strats = [TARGETS_BY_STRATEGY[0]]
    else:
        strats = [TARGETS_BY_STRATEGY[1]]
    return map(lambda x: makeTarget(*x), strats)
    
targets = []

def findTargets(names):
    def find1(name):
        for target in targets:
            if target.name == name:
                return target
        return None
    return map(find1, names)

def matches(string, pattern):
    return string.find(pattern) != -1

def main():
    global OPTIONS, targets
    usage = "usage: %prog stage [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("", "--dry", dest="dry",
                        action='store_true', default=False,
                        help="If provided, nothing actually gets run, just a dry run")
    parser.add_option("-v", "--verbose", dest="verbose",
                        action='store_true', default=False,
                        help="If provided, print out a lot of information")
    parser.add_option("", "--byQEval", dest="byQEval",
                        action='store_true', default=False,
                        help="If provided, variant eval will be run by Q threshold")
    parser.add_option("-s", "--splitByChr", dest="splitByChr",
                        action='store_true', default=False,
                        help="If provided, we'll parallelize by chromosome over the farm")
    parser.add_option("", "--dev", dest="dev",
                        action='store_true', default=False,
                        help="If provided, we'll use the GATK dev build")
    parser.add_option("-d", "--dir", dest="dir",
                        type='string', default="",
                        help="If provided, this is the root where files are read and written")
    parser.add_option("-L", "", dest="WGSIntervals",
                        type='string', default=None,
                        help="If provided, these are the interval files we will process for WGS")
    parser.add_option("-q", "--farm", dest="farmQueue",
                        type="string", default=None,
                        help="Farm queue to send processing jobs to")
    parser.add_option("-p", "--parallel", dest="parallel",
                        type="int", default=None,
                        help="Number of parallel shared memory threads")
    parser.add_option("-t", "--target", dest="target",
                        type="string", default=None,
                        help="Only run jobs with names containing this string")
    parser.add_option("", "--noRaw", dest="noRaw",
                        action="store_true", default=False,
                        help="Exclude raw calls from output")
                       
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    # set up targets
    # WGS
    WGS_INTERVAL = OPTIONS.WGSIntervals
    #targets += targetsByStrategy(CallTarget, 'GA2.WGS.cleaned', BAM_ROOT + '/NA12878.GA2.WGS.bwa.cleaned.bam', WGS_INTERVAL, filters = WGS_FILTER)
    #targets.append(CallTarget('GA2.WGS.raw', '/seq/dirseq/pem/seq/picard_aggregation/G2946gaII/NA12878/v1/NA12878.bam', WGS_INTERVAL, filters = WGS_FILTER))

    # HiSeq
    targets += targetsByStrategy(CallTarget, 'HiSeq.WGS.raw', '/seq/dirseq/pem/seq/picard_aggregation/G2946/NA12878/v1/NA12878.bam', WGS_INTERVAL, filters = WGS_FILTER)
    #targets += targetsByStrategy(CallTarget, 'HiSeq.WGS.cleaned', '/humgen/1kg/analysis/bamsForDataProcessingPapers/scriptsToMakeBams/tmp.list', WGS_INTERVAL, filters = WGS_FILTER)
    targets += targetsByStrategy(CallTarget, 'HiSeq.WGS.cleaned', BAM_ROOT + '/NA12878.HiSeq.WGS.bwa.cleaned.recal.bam', WGS_INTERVAL, filters = WGS_FILTER)
    
    # WE
    targets += targetsByStrategy(weTarget, 'GA2.WEx.cleaned', BAM_ROOT + '/NA12878.WEx.cleaned.recal.bam')
    targets += targetsByStrategy(weTarget, 'GA2.WEx.raw', '/seq/picard_aggregation/C308/NA12878/v3/NA12878.bam')
    #targets.append(weTarget('GA2.WEx.raw', '/seq/picard_aggregation/C308/NA12878/v3/NA12878.bam'))

    #targets += targetsByStrategy(CallTarget, 'CG.WGS.raw', '/seq/complete_genomics/GS00106-DNA_E01-180_NA12878/SAM0/merge/NA12878.bam', WGS_INTERVAL, filters = WGS_FILTER)

    # CG
    # todo -- fixme -- needs genome-wide bams on hg18
    #targets.append(CallTarget('CG.chr1.raw', '/humgen/gsa-hphome1/kiran/one_off_projects/multiTechComparisons/results/v5/NA12878/CG.full/sample.chr1.primaryAlignmentsMarked.dupesRemoved.bam', WGS_INTERVAL, filters = WGS_FILTER, callArgs = "-bm THREE_STATE", b36 = True))

    # low-pass
    # '/humgen/gsa-hpprojects/1kg/1kg_pilot1/freeze5_merged/low_coverage_CEU.1.bam
    # targets.append(CallTarget('CEU.lowpass.cleaned', 'CEU.bam.list', WGS_INTERVAL, b36 = True)) 

    # 1KG SLX
    # targets.append(CallTarget('1KG.NA12878', '/humgen/gsa-hpprojects/1kg/1kg_pilot2/useTheseBamsForAnalyses/NA12878.SLX.bam', WGS_INTERVAL, b36 = True)) 

    # MCDK1 special case
    #MCKD1_INTERVAL = "chr1:152448527-154998173"
    #targets.append(CallTarget('MCKD1.raw', '/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/MCKD1_WGS/MCKD1.bam.list', MCKD1_INTERVAL, filters = WGS_FILTER, callArgs = '-bm THREE_STATE')) 
    #targets.append(CallTarget('MCKD1.cleaned', '/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/MCKD1_WGS/bams_linkage/MCKD1.bam.cleaned.bam.chr1:152448527-154998173.bam', MCKD1_INTERVAL, filters = WGS_FILTER, callArgs = '-bm THREE_STATE')) 

    stages = map(string.lower, args[0].split(","))
    STAGES = ['callsnps', 'callindels', 'indelmask', 'snpfilter', 'indelfilter', 'to_hg18', 'optimize', 'eval', 'confusion_matrix']
    for stage in stages:
        if stage not in STAGES:
            sys.exit('unknown stage ' + stage)

    if OPTIONS.dir != "" and not os.path.exists(OPTIONS.dir):
        os.makedirs(OPTIONS.dir)

    allJobs = []
    
    def includeStage(name):
        return name in stages

    for callTarget in targets:
        if "raw" in callTarget.name and OPTIONS.noRaw:
            print 'Skipping raw data', callTarget
            continue

        # setup pipeline args
        GATK_JAR = GATK_STABLE_JAR
        if ( OPTIONS.dev ): GATK_JAR = GATK_DEV_JAR
        myPipelineArgs = PipelineArgs(GATK_JAR = GATK_JAR, name = callTarget.name, excludeChrs = EXCLUDE_CHRS)
        myPipelineArgs.addGATKArg(EXTRA_GATK_ARGS)
        myPipelineArgs.addGATKArg(callTarget.getIntervalArg())
        if ( OPTIONS.parallel != None ): myPipelineArgs.addGATKArg(' -nt ' + OPTIONS.parallel)
        if ( callTarget.b36 ): myPipelineArgs.ref = 'b36'
        print callTarget.name, callTarget.b36, myPipelineArgs.ref
    
        lastJobs = None
        if OPTIONS.target != None and not matches(callTarget.name, OPTIONS.target):
            print 'Skipping target', callTarget
            continue
        
        def updateNewJobs(newjobs, lastJobs):
            if OPTIONS.verbose:
                print 'New jobs', newjobs
            #for job in newjobs:
            #    print '    job ', job
            allJobs.append(newjobs)
            if newjobs != []: 
                lastJobs = newjobs
            return [], lastJobs

        newJobs = []
        def execStage(name, func, vcf = None, lastJobs = []):
            if OPTIONS.verbose: print 'Name is', name
            newJobs, newVcf = func(myPipelineArgs, callTarget, vcf, lastJobs)
            if newVcf != None: vcf = newVcf
            if OPTIONS.verbose: print 'VCF is', vcf
            callTarget.addVcf(vcf)
            if includeStage(name): newJobs, lastJobs = updateNewJobs(newJobs, lastJobs)
            if OPTIONS.verbose: print 'execStage:', newJobs, lastJobs, vcf
            return newJobs, lastJobs, vcf

            
        newJobs, callSNPJobs, vcf = execStage('callsnps', callSNPs)
        newJobs, lastJobs, vcf = execStage('to_hg18', convertToHg18, vcf, callSNPJobs)
        newJobs, filterSNPJobs, vcf = execStage('snpfilter', filterSNPs, vcf, callSNPJobs)

        # indel jobs
        newJobs, callIndelJobs, vcf = execStage('callindels', callIndels, vcf)
        newJobs, indelMaskJobs, vcf = execStage('indelmask', createIndelMask, vcf, callIndelJobs)
        newJobs, filterIndelJobs, vcf = execStage('indelfilter', filterIndels, vcf, indelMaskJobs + filterSNPJobs)

        # optimization
        newJobs, optimizeJobs, vcf = execStage('optimize', VariantOptimizer, vcf, filterIndelJobs)

        # eval
        newJobs, evalJobs, vcf = execStage('eval', evalSNPs, vcf, optimizeJobs)

        # newJobs, lastJobs, ignore = execStage('confusion_matrix', computeConfusionMatrix, vcf)

    print 'EXECUTING JOBS'
    executeJobs(allJobs, farm_queue = OPTIONS.farmQueue, just_print_commands = OPTIONS.dry) 

#
# Actual commands
#
def convertToHg18( myPipelineArgs, callTarget, vcf, lastJobs ):
    if callTarget.b36:
        outputVCF = vcf.replace(".b36", "")
        cmd = 'python /humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/python/vcf_b36_to_hg18.py ' + vcf + ' ' + outputVCF  
        jobs = [FarmJob(cmd, jobName = callTarget.name + '.' + 'b36ToHg18', dependencies = lastJobs)]
        return jobs, outputVCF
    else:
        return [], vcf

def callSNPs( myPipelineArgs, callTarget, ignore, lastJobs ):
    outputVCF = appendExtension(callTarget.getVcf(), "ug")
    if callTarget.b36:
        outputVCF = appendExtension(outputVCF, "b36")
    print 'outputVCF', outputVCF
    ugArgs = '-T UnifiedGenotyperV2 -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_hg18.rod -I %s %s -o %s %s' % (callTarget.getBam(), UG_ARGS, outputVCF, callTarget.getCallArgs())
    farmCmds = simpleGATKCommand( myPipelineArgs, 'call.UG', ugArgs, lastJobs )
    if OPTIONS.splitByChr and not callTarget.hasInterval():
        farmCmds = splitGATKCommandByChr( myPipelineArgs, farmCmds[0], [outputVCF], [mergeVCFs] )
    return farmCmds, outputVCF

INDEL_MASK_SIZE = 10
def getIndelCallFiles(callTarget):
    outputBed = appendExtension(callTarget.getVcf(), "indels.bed", False)
    outputVerbose = appendExtension(callTarget.getVcf(), "indels.verbose.txt", False)
    outputVCF = appendExtension(callTarget.getVcf(), "indels.vcf", False)
    outputMask = appendExtension(callTarget.getVcf(), "indels.%d.mask" % INDEL_MASK_SIZE, False)
    return outputBed, outputVerbose, outputVCF, outputMask

def callIndels( myPipelineArgs, callTarget, ignore, lastJobs ):
    outputBed, outputVerbose, indelsVCF, outputMask = getIndelCallFiles(callTarget)
    IGV2_ARGS = '-T IndelGenotyperV2 -ws 500 -I %s -bed %s -verbose %s -o %s -rf Platform454' % (callTarget.getBam(), outputBed, outputVerbose, indelsVCF)
    farmCmds = simpleGATKCommand( myPipelineArgs, 'call.CallIndels', IGV2_ARGS, lastJobs )
    if OPTIONS.splitByChr and not callTarget.hasInterval():
        farmCmds = splitGATKCommandByChr( myPipelineArgs, farmCmds[0], [outputBed, outputVerbose], [mergeByCat, mergeByCat] )
    return farmCmds, None

def createIndelMask( myPipelineArgs, callTarget, vcf, lastJobs ):
    outputBed, outputVerbose, outputVCF, outputMask = getIndelCallFiles(callTarget)
    cmd = 'python /humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/python/makeIndelMask.py %s %d %s' % (outputBed, INDEL_MASK_SIZE, outputMask)  
    jobs = [FarmJob(cmd, jobName = callTarget.name + '.' + 'makeIndelMask', dependencies = lastJobs)]
    return jobs, None

def filterSNPs(myPipelineArgs, callTarget, vcf, lastJobs ):
    out = appendExtension(vcf, 'snpfiltered')
    filterString = ' '.join(map(lambda x: '--filterName %s --filterExpression "%s"' % (x[0], x[1]), callTarget.filters))
    return simpleGATKCommand( myPipelineArgs, 'call.filterSNPs', '-T VariantFiltration -B:variant,VCF %s -o %s --filterName LowQual --filterExpression "QUAL < 50.0" --clusterWindowSize 10 --filterName HARD_TO_VALIDATE --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" %s' % ( vcf, out, filterString ), lastJobs ), out

def filterIndels(myPipelineArgs, callTarget, vcf, lastJobs ):
    out = appendExtension(vcf, 'indelfiltered')
    outputBed, outputVerbose, outputVCF, outputMask = getIndelCallFiles(callTarget)
    return simpleGATKCommand( myPipelineArgs, 'call.filterIndels', '-T VariantFiltration -B:variant,VCF %s -o %s --maskName Indel -B:mask,Bed %s' % ( vcf, out, outputMask ), lastJobs ), out

def VariantOptimizer( myPipelineArgs, callTarget, vcf, lastJobs ):
    if callTarget.optimize:
        clusterFile = appendExtension(vcf, 'optimized.clusters', False)
        tranchesFile = appendExtension(vcf, 'optimized.tranches', False)
        optOutVCF = appendExtension(vcf, 'optimized')
        out = appendExtension(vcf, 'optimized.cut')
        table = appendExtension(vcf, 'optimized.table', False)

        #NCLUSTERS = 4
        #ITERATIONS = 3
        #ITERATION_TO_TAKE = ITERATIONS
        DBSNP_PRIOR = 2.0 # dbSNP seems dodger and dodger        
        IGNORE_FILTERS_CLUSTERING = "-ignoreFilter DPFilter -ignoreFilter ABFilter -ignoreFilter ESPStandard"
        #IGNORE_FILTERS_CLUSTERING = "-ignoreFilter DPFilter -ignoreFilter ABFilter -ignoreFilter LowQual -ignoreFilter ESPStandard"
        IGNORE_FILTERS_SCORING = IGNORE_FILTERS_CLUSTERING + " -ignoreFilter HARD_TO_VALIDATE"
        annotationsToOptimize = ['SB', 'HaplotypeScore', "QD", 'HRun']
        annotationsToOptimizeArg = ' '.join(map(lambda x: '-an ' + x, annotationsToOptimize)) # '' ['DP', 'SB', 'HaplotypeScore', 'MQ', "QD", 'HRun']
        #tranches = ' '.join(map( lambda x: '-tranche ' + str(x), [1, 5, 10]))
        tranches = ' '.join(map( lambda x: '-tranche ' + str(x), [0.1, 1, 10]))
        maxVariantsToShow = 2500
        #singletonFPRate = 0.2
        #hapmapVCF = '/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.2/genotypes_r27_nr.hg18_fwd.vcf'
        hapmapVCF = 'hapmap_analysis/sitesr27_nr.hg18_fwd.vcf'

        REGENERATE_VARIANT_CLUSTERS = True
        if ( REGENERATE_VARIANT_CLUSTERS ):
            jobs1 = simpleGATKCommand( myPipelineArgs, 'call.GenerateVariantClusters', '-T GenerateVariantClusters -qual %d -std 3.5 -mG %d -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_hg18.rod -B:hapmap,vcf %s -B:input,VCF %s -clusterFile %s %s %s --NoByHapMapValidationStatus' % ( callTarget.minQual, callTarget.maxClusters, hapmapVCF, vcf, clusterFile, annotationsToOptimizeArg, IGNORE_FILTERS_CLUSTERING ), lastJobs )
            #jobs1 = simpleGATKCommand( myPipelineArgs, 'call.GenerateVariantClusters', '-T GenerateVariantClusters -qual %d -std 3.5 -mG %d -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_hg18.rod -B:input,VCF %s -clusterFile %s %s %s' % ( callTarget.minQual, callTarget.maxClusters, vcf, clusterFile, annotationsToOptimizeArg, IGNORE_FILTERS_CLUSTERING ), lastJobs )
        else:
            jobs1 = lastJobs
        jobs2 = simpleGATKCommand( myPipelineArgs, 'call.VariantRecalibrator', '-T VariantRecalibrator -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_hg18.rod -B:input,VCF %s -clusterFile %s -o %s --target_titv %f %s -resources ~/dev/GenomeAnalysisTK/trunk/R/ %s -priorDBSNP %.2f -tranchesFile %s' % ( vcf, clusterFile, optOutVCF, callTarget.targetTiTv, IGNORE_FILTERS_SCORING, tranches, DBSNP_PRIOR, tranchesFile ), jobs1 )

        cmd21 = 'python /humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/python/vcf2table.py -f CHROM,POS,ID,AC,AF,AN,DB,' + ','.join(annotationsToOptimize) + ' ' + vcf + ' -o ' + table
        jobs21 = [FarmJob(cmd21, jobName = callTarget.name + '.call.' + 'VariantRecalibrationReport.vcf2table', dependencies = jobs2)]

        cmd22 = 'Rscript /humgen/gsa-scr1/depristo/dev/GenomeAnalysisTK/trunk/R/VariantRecalibratorReport/VariantRecalibratorReport.R %s %s %s NA %d' % (clusterFile, clusterFile, table, maxVariantsToShow)
        jobs22 = [FarmJob(cmd22, jobName = callTarget.name + '.call.' + 'VariantRecalibrationReport.RScript', dependencies = jobs21)]

        jobs3 = simpleGATKCommand( myPipelineArgs, 'call.ApplyVariantCuts', '-T ApplyVariantCuts -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_hg18.rod -B:input,VCF %s -o %s -tranchesFile %s --fdr_filter_level %f' % ( optOutVCF, out, tranchesFile, callTarget.tranchToTake ), jobs2 )
        return jobs1 + jobs2 + jobs21 + jobs22 + jobs3, out
    else:
        return [], vcf

def computeConfusionMatrix(myPipelineArgs, callTarget, vcf, lastJobs ):
    out = appendExtension(vcf, 'confusionmatrix', addExtension=False)
    CM_ARGS = '-T ComputeConfusionMatrix -I %s -o %s' % (callTarget.getBam(), out)
    farmCmds = simpleGATKCommand( myPipelineArgs, 'ConfusionMatrix', CM_ARGS, lastJobs )
    return farmCmds, None

def evalSNPs(myPipelineArgs, callTarget, vcf, lastJobs):
    evalRoot = OPTIONS.dir
 
    oldMemory = myPipelineArgs.memory
    myPipelineArgs.memory = '3g'
    def eval1(vcf, namePostfix = "", args = ""):
        out = os.path.join(OPTIONS.dir, os.path.basename(vcf) + namePostfix + ".ve2")
        maybeHiSeqBindings = ""
        hiSeqComp = os.path.join(OPTIONS.dir,"HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.optimized.cut.vcf")
        omni = " -B:compOmni,VCF Omni.NA12878.hg18.vcf"
        if os.path.exists(hiSeqComp):
            maybeHiSeqBindings = "-B:comp_HiSeq,VCF " + hiSeqComp + " "
        validation_bindings = maybeHiSeqBindings + "-B:comp_p2_val,VCF 1kg_pilot2_snps.hg18.vcf -B:comp_CG,VCF /humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/NA12878/CG.hg18.vcf -B:compTrio,VCF CEU.trio.2010_03.genotypes.vcf -B:compTrioNovel,VCF CEU.trio.novels.2010_03.genotypes.vcf -B:comp_hm3,VCF /humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.2/by_population/genotypes_CEU_phase3.2_consensus.hg18_fwd.vcf " + omni # not in hg18 space :-(
        tranches = ""
        if vcf.find("optimized") != -1:
            args += " -tf " + appendExtension(vcf.replace(".cut", ""), 'tranches', False)
            vcf = appendExtension(vcf.replace(".cut", ""), 'vcf', False)
        gatk_args = ("-T VariantEval -reportType Grep -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_hg18.rod -B:eval,VCF %s " + validation_bindings + " -sample NA12878 -o %s -E CompOverlap -E GenotypeConcordance -E TiTvVariantEvaluator -E CountVariants %s") % ( vcf, out, args )

        name = "EVAL_%s_%s" % (callTarget.name, namePostfix)
        return simpleGATKCommand( myPipelineArgs, name, gatk_args, lastJobs )[0]
 
    jobs = []
    for vcf in callTarget.getVcfs():
        jobs.append(eval1(vcf))
        if OPTIONS.byQEval and vcf.find("optimized") != -1:
            for Q in [0.01, 0.02, 0.03, 0.03]:
                jobs.append(eval1(vcf, '.Q' + str(Q), '-Q ' + str(Q)))

    myPipelineArgs.memory = oldMemory
    return jobs, None

if __name__ == "__main__":
    main()
