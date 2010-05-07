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

def main():
    global OPTIONS
    usage = "usage: %prog [options] stages input.bam outputRoot"
    parser = OptionParser(usage=usage)
    parser.add_option("", "--dry", dest="dry",
                        action='store_true', default=False,
                        help="If provided, nothing actually gets run, just a dry run")
    parser.add_option("", "--b36", dest="useB36",
                        action='store_true', default=False,
                        help="If provided, BAM is assumed to aligned to b36 named chromosomes")
    parser.add_option("-v", "--verbose", dest="verbose",
                        action='store_true', default=False,
                        help="")
    parser.add_option("-n", "--name", dest="name",
                        type="string", default="realignBamByChr",
                        help="Farm queue to send processing jobs to")
    parser.add_option("-q", "--farm", dest="farmQueue",
                        type="string", default=None,
                        help="Farm queue to send processing jobs to")
    parser.add_option("-e", "--extraArgs", dest="extraArgs",
                        type="string", default=None,
                        help="")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")

    stages = map(string.lower, args[0].split(","))
    inputBam, outputRoot = args[1:]
    outputBamList = outputRoot + '.bams.list'
    
    STAGES = ['targets', 'realign', 'index']
    for stage in stages:
        if stage not in STAGES:
            sys.exit('unknown stage ' + stage)

    myPipelineArgs = PipelineArgs(name = OPTIONS.name)
    if ( OPTIONS.useB36 ):
        myPipelineArgs.ref = 'b36'

    allJobs = []
    
    def includeStage(name):
        return name in stages

    out = open(outputBamList, 'w')
    for chr in hg18:
        lastJobs = None
        
        def updateNewJobs(newjobs, lastJobs):
            if OPTIONS.verbose:
                print 'New jobs', newjobs
            allJobs.append(newjobs)
            if newjobs != []: 
                lastJobs = newjobs
            return [], lastJobs

        newJobs = []

        def execStage(name, func, args = [], lastJobs = []):
            if OPTIONS.verbose: print 'Name is', name
            newJobs, results = func(myPipelineArgs, chr, inputBam, outputRoot + '.' + chr, args, lastJobs)
            if includeStage(name): newJobs, lastJobs = updateNewJobs(newJobs, lastJobs)
            return newJobs, lastJobs, results
            
        newJobs, lastJobs, intervals = execStage('targets', createTargets)
        newJobs, lastJobs, realignedBam = execStage('realign', realign, intervals, lastJobs)
        # need to merge and then index
        newJobs, lastJobs, ignore = execStage('index', index, realignedBam, lastJobs)
        print >> out, os.path.abspath(realignedBam)
     
    out.close() 
    print 'EXECUTING JOBS'
    executeJobs(allJobs, farm_queue = OPTIONS.farmQueue, just_print_commands = OPTIONS.dry) 

def createTargets( myPipelineArgs, chr, inputBam, outputRoot, args, lastJobs ):
    outputIntervals = outputRoot + ".intervals"
    GATKArgs = '-T RealignerTargetCreator -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_hg18.rod -I %s -o %s -mrl 10000 -L %s' % (inputBam, outputIntervals, chr)
    return simpleGATKCommand( myPipelineArgs, 'CreateInterval' + chr, GATKArgs, lastJobs ), outputIntervals

def realign( myPipelineArgs, chr, inputBam, outputRoot, intervals, lastJobs ):
    outputBAM = outputRoot + ".bam"
    GATKArgs = '-T IndelRealigner -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_hg18.rod -I %s -targetIntervals %s --output %s -sort ON_DISK -mrl 100000 -L %s' % (inputBam, intervals, outputBAM, chr)
    return simpleGATKCommand( myPipelineArgs, 'Realign' + chr, GATKArgs, lastJobs ), outputBAM

def index( myPipelineArgs, chr, inputBam, outputRoot, realignedBam, lastJobs ):
    return indexBAMFile( myPipelineArgs.name, realignedBam, lastJobs )

if __name__ == "__main__":
    main()
