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
    usage = "usage: %prog [options] indels.verbose.txt"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", dest="verbose",
                        action='store_true', default=False,
                        help="")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    indelsFile = args[0]

    print 'event size nReadsWithIndel nReadsWithoutIndel nReadsTotal'
    for line in open(indelsFile):
        try:
            parts = line.split()
            # chr1    30502   30505   -TTT    OBS_COUNTS[C/A/T]:5/5/13        AV_MM[C/R]:0.00/0.75    AV_MAPQ[C/R]:37.00/27.13        NQS_MM_RATE[C/R]:0.00/0.01      NQS_AV_QUAL[C/R]:30.90/28.30    STRAND_COUNTS[C/C/R/R]:2/3/6/2
            chr, start, stop = parts[0:3]
            event = parts[3]
            counts = parts[4]
            isInsertion = event[0] == '+'
            size = len(event[1:])
            nConsensusReads, nReadsWithIndel, nReads = map(int, counts.split(':')[1].split('/'))
            nReadsWithoutIndel = nReads - nReadsWithIndel
            print event[0], size, nReadsWithIndel, nReadsWithoutIndel, nReads     
        except:
            continue
            
if __name__ == "__main__":
    main()
