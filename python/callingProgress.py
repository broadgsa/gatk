import farm_commands
import os.path
import sys
from optparse import OptionParser
from datetime import date
import glob
import operator

if __name__ == "__main__":
    usage = "usage: %prog files"
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--size", dest="regionSize",
                        type="int", default=None, help="")
                         
    (OPTIONS, args) = parser.parse_args()

    def extract(line):
        s = line.split()
        t = s[0].split(":")
        if len(t) == 2:
            chr, pos = t
            return chr, int(pos)
        else:
            return None, None

    for file in args:
        chr, lastPos, startPos, calledBases, progress = None, None, None, None, None
        lastLine = None
        for line in open(file):
            if startPos == None: 
                chr, startPos = extract(line)
            lastLine = line
        if lastLine <> None and startPos <> None:
            lastPos = extract(lastLine)[1]
            calledBases = lastPos - startPos
            progress = "%.2f" % (float(calledBases) / OPTIONS.regionSize * 100)
        print file, chr, startPos, lastPos, calledBases, progress