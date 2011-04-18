#!/bin/tcsh
import os.path
import sys
from optparse import OptionParser
from itertools import *

def main():
    global OPTIONS
    usage = "usage: %prog [options] mode hg19Tohg18.chain"
    parser = OptionParser(usage=usage)
    
#     parser.add_option("-D", "--delete_while_archiving", dest="reallyDeleteInArchiveMode",
#                         action='store_true', default=False,
#                         help="if provided, we'll actually delete records when running in archive mode")
         
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("Requires exact 1 chain to analyze")

    hg192hg18 = args[0]
    
    writeChain(hg192hg18, "b37tohg18.chain", lambda x: hg2b(x, 2))
    writeChain(hg192hg18, "b37tob36.chain", lambda x: hg2b(hg2b(x, 2), 7))
    
HG2BCONTIG = dict()
for c in range(1, 23) + ["X", "Y"]:
    HG2BCONTIG["chr" + str(c)] = str(c)
HG2BCONTIG["chrM"] = "MT"
    
def hg2b(line, pos):
    parts = line.split()
    if len(parts) > pos and "chr" in parts[pos]:
        if parts[pos] in HG2BCONTIG:
            parts[pos] = HG2BCONTIG[parts[pos]]
        else:
            print 'Skipping ', parts[pos]
    return "\t".join(parts)
    
def writeChain(inFile, outFile, transform):
    out = open(outFile, "w")
    for line in open(inFile):
        newLine = transform(line)
        #print 'newline', newLine
        if newLine != None:
            out.write(newLine)
            out.write("\n")
    out.close()

if __name__ == "__main__":
    main()
