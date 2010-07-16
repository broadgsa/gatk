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

GATK_STABLE_JAR = '/home/radon01/depristo/dev/GenomeAnalysisTKStable/trunk/dist/GenomeAnalysisTK.jar'
GATK_DEV_JAR = '/home/radon01/depristo/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar'
GATK_JAR = GATK_STABLE_JAR

# add to GATK to enable dbSNP aware cleaning
# -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_hg18.rod

hg18 = ['chrM'] + ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']
b36 = [str(i) for i in range(1,23)] + ['X', 'Y', 'MT']

#hg18 = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']
#b36 = [str(i) for i in range(1,23)] + ['X', 'Y']

HG18_TO_B36 = {
    'hg18' : 'b36',
    '/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta' : '/broad/1KG/reference/human_b36_both.fasta',
    '/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta.fai' : '/broad/1KG/reference/human_b36_both.fasta.fai',
    'chrM' : 'MT',
    'chr' : '' }
    
def hg18args_to_b36(s):
    for key, value in HG18_TO_B36.iteritems():
         s = s.replace(key, value)
    return s

def appendExtension(path, newExt, addExtension = True):
    root, basename = os.path.split(path)
    basename, ext = os.path.splitext(basename)
    print root, basename, ext, path
    s = basename + '.' + newExt
    if addExtension: s += ext
    return os.path.join(root, s)
#    return os.path.join(OPTIONS.dir, s)

class PipelineArgs:
    def __init__( self, GATK_JAR = GATK_JAR, ref = 'hg18', name = None, memory = '3g', excludeChrs = [] ):
        self.GATK = 'java -Xmx%s -Djava.io.tmpdir=/broad/shptmp/depristo/tmp/ -jar ' + GATK_JAR + ' -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -l INFO '
        self.ref = ref
        self.name = name
        self.memory = memory
        self.excludeChrs = excludeChrs

    def convertToB36(self): 
        return self.ref == 'b36'

    def addGATKArg(self, arg):
        if len(arg) > 0 and arg[0] != ' ':
            arg = ' ' + arg
        self.GATK += arg

    def getCommandName(self, suffix):
        if self.name != None:
            return self.name + "." + suffix
        else:
            return suffix
            
    def finalizedGATKCommand(self, args):
        cmd = (self.GATK % self.memory) + ' ' + args
        if self.convertToB36():
            cmd = hg18args_to_b36(cmd)
        return cmd
        
    def chrsToSplitBy(self, chrs):
        #print 'XXX', self.excludeChrs
        return filter(lambda x: x not in self.excludeChrs, chrs)
        
# 
# General features
#
def simpleGATKCommand( pargs, name, args, lastJobs ):
    cmd = pargs.finalizedGATKCommand(args)
    return [FarmJob(cmd, jobName = pargs.getCommandName(name), dependencies = lastJobs)]

# 
# Takes a simpleGATKCommand and splits it by chromosome, merging output
# 
def splitGATKCommandByChr( myPipelineArgs, cmd, outputsToParallelize, mergeCommands ):
    if cmd.cmd_str_from_user.find(" -L") != -1:
        sys.exit("Found -L argument in command -- cannot be provided for parallelization by chromosome: " + cmd.cmd_str_from_user)
    def makeChrCmd(chr):
        chrOutputMap = map(lambda x: appendExtension(x, chr), outputsToParallelize)
        chr_cmd_str = cmd.cmd_str_from_user 
        chr_cmd_str += ' -L ' + chr
        for x, y in zip(outputsToParallelize, chrOutputMap):
            chr_cmd_str = chr_cmd_str.replace(x, y)

        if myPipelineArgs.convertToB36():
            chr_cmd_str = hg18args_to_b36(chr_cmd_str)
            
        chrCmd = FarmJob(chr_cmd_str, jobName = cmd.jobName + '.byChr' + chr, dependencies = cmd.dependencies)
        return chrCmd, chrOutputMap

    #print '######################################### chrsToSplitBy', myPipelineArgs.chrsToSplitBy(hg18)
    splits = map( makeChrCmd, myPipelineArgs.chrsToSplitBy(hg18) )
    splitCommands = map(lambda x: x[0], splits)

    def mergeCommand1(i):
        mergeCommand = mergeCommands[i]
        mergeFile = outputsToParallelize[i]
        splitFiles = map(lambda x: x[1][i], splits)
        return FarmJob(mergeCommand(splitFiles, mergeFile), jobName = cmd.jobName + '.merge', dependencies = splitCommands)
        
    mergeCommands = map(mergeCommand1, range(len(outputsToParallelize)))

    return splitCommands + mergeCommands

def mergeVCFs(splitFiles, mergeFile):
    splitFilesString = ' '.join(splitFiles)
    cmd = "python ~/dev/GenomeAnalysisTK/trunk/python/mergeVCFs.py -a -f /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta.fai %s > %s" % (splitFilesString, mergeFile)
    return cmd

def mergeByCat(splitFiles, mergeFile):
    splitFilesString = ' '.join(splitFiles)
    return "cat %s > %s" % (splitFilesString, mergeFile)

def indexBAMFile( name, bamFile, lastJobs ):
    cmd = 'samtools index ' + bamFile
    return [FarmJob(cmd, jobName = 'samtools.index.' + name, dependencies = lastJobs)], None
