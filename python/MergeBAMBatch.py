import farm_commands
import os.path
import sys
from optparse import OptionParser
from datetime import date
import glob
import operator
import ValidateGATK
import picard_utils
from MergeBAMsUtils import *

def splitSourcesByPopulation(allSources, merged_filename_base, NAID2Pop):
    sourcePairs = [[source, source] for source in allSources]
    return groupSources(sourcePairs, NAID2Pop, merged_filename_base)

if __name__ == "__main__":
    usage = "usage: %prog files.list [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-q", "--farm", dest="farmQueue",
                        type="string", default=None,
                        help="Farm queue to send processing jobs to")
    parser.add_option("-d", "--dir", dest="output_dir",
                        type="string", default="./",
                        help="Output directory")
    parser.add_option("", "--dry", dest="dry",
                        action='store_true', default=False,
                        help="If provided, nothing actually gets run, just a dry run")
    parser.add_option("-i", "--ignoreExistingFiles", dest="ignoreExistingFiles",
                        action='store_true', default=False,
                        help="Ignores already written files, if present")
    parser.add_option("-s", "--useSamtools", dest="useSamtools",
                        action='store_true', default=False,
                        help="If present, uses samtools to perform the merge")
    parser.add_option("-m", "--mergeBin", dest="mergeBin",
                        type="string", default=None,
                        help="Path to merge binary")
    parser.add_option("-n", "--naIDPops", dest="NAIDS2POP",
                        type="string", default=None,
                        help="Path to file contains NAID POP names.  If provided, input files will be merged by population")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    directory = OPTIONS.output_dir
    
    if not os.path.exists(directory):
        os.mkdir(directory)
    
    NAID2Pop = None
    if OPTIONS.NAIDS2POP <> None:
        NAID2Pop = readNAIdMap(OPTIONS.NAIDS2POP)
    
    for line in open(args[0]):
        s = line.split()
        if ( s <> [] and s[0] <> '#' ):
            merged_filename_base = s[0]
            allSources = reduce( operator.__add__, map( glob.glob, s[1:] ), [] )
            print 'Merging info:'
            for spec in splitSourcesByPopulation(allSources, merged_filename_base, NAID2Pop):
                spec.setPath(directory)
                spec.pprint()
                
                jobid = None
                if OPTIONS.ignoreExistingFiles or not os.path.exists(spec.getMergedBAM()):
                    output = spec.getMergedBase() + '.stdout'
                    cmd = spec.mergeCmd(OPTIONS.mergeBin, useSamtools = OPTIONS.useSamtools)
                    #print cmd
                    jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, output, just_print_commands = OPTIONS.dry)
    
                if OPTIONS.ignoreExistingFiles or not os.path.exists(spec.getMergedBAMIndex()):
                    jobid = farm_commands.cmd(spec.getIndexCmd(), OPTIONS.farmQueue, '', waitID = jobid, just_print_commands = OPTIONS.dry)

