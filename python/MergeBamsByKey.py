import farm_commands
import os.path
import sys
from optparse import OptionParser
import picard_utils
from MergeBAMsUtils import *

def splitSourcesByKeys( bams, keys ):
    keyPairs = [[key, key] for key in keys]
    keybamPairs = zip(keys, bams)
    return groupSources(keybamPairs, keyPairs, None)

if __name__ == "__main__":
    usage = """usage: %prog bams.list [options]
Merges BAM files by keys from a file of a list of bams.

bams.list is a whitespace separated file.  One column (--keyCol arg) is the key, and another
column (--bamCol) is a path to a bam file.  This program will group the bam files
by key and spawn merge and index jobs to merge all of the files sharing the same key together"""

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
    parser.add_option("-m", "--mergeBin", dest="mergeBin",
                        type="string", default=None,
                        help="Path to merge binary")
    parser.add_option("", "--keyCol", dest="keyCol",
                        type=int, default=1,
                        help="Column in the list file holding the key")
    parser.add_option("", "--bamCol", dest="bamCol",
                        type=int, default=2,
                        help="Column in the list file holding the bam file path")
    parser.add_option("-l", "--link", dest="link",
                        action='store_true', default=False,
                        help="If true, program will soft link single bam files that don't need merging")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    directory = OPTIONS.output_dir
    
    if not os.path.exists(directory):
        os.mkdir(directory)
    
    bamsList = [line.strip().split() for line in open(args[0])]
    keys = map( lambda x: x[OPTIONS.keyCol-1], bamsList ) 
    bams = map( lambda x: x[OPTIONS.bamCol-1], bamsList ) 
    
    print 'Merging info:'
    for info in bamsList: print info
    for spec in splitSourcesByKeys(bams, keys):
        spec.setPath(directory)
        spec.pprint()
        
        jobid = None
        if OPTIONS.ignoreExistingFiles or not os.path.exists(spec.getMergedBAM()):
            output = spec.getMergedBase() + '.stdout'
            if len(spec.sources()) == 1 and OPTIONS.link:
                cmd = 'ln -s ' + spec.sources()[0] + ' ' + spec.getMergedBAM()
            else:
                cmd = spec.mergeCmd(OPTIONS.mergeBin)
            print cmd
            jobid = farm_commands.cmd(cmd, OPTIONS.farmQueue, output, just_print_commands = OPTIONS.dry)

        if OPTIONS.ignoreExistingFiles or not os.path.exists(spec.getMergedBAMIndex()):
            pass
            jobid = farm_commands.cmd(spec.getIndexCmd(), OPTIONS.farmQueue, '', waitID = jobid, just_print_commands = OPTIONS.dry)

