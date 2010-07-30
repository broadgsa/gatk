from farm_commands2 import *
import os
import sys
from optparse import OptionParser
from datetime import date
import glob
import operator
import faiReader
import math
import shutil
import string

def main():
    global OPTIONS
    usage = "usage: %prog [options] b36VCF hg18VCF"
    parser = OptionParser(usage=usage)
    parser.add_option("", "--dry", dest="dry",
                        action='store_true', default=False,
                        help="If provided, nothing actually gets run, just a dry run")
                       
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")

    b36vcf, hg18vcf = args

    temp = open("tmp", 'w')
    mitotemp = open("mtmp",'w')
    for line in open(b36vcf):
        length = len(line)
        if length > 2 and line[0] != '#':
            if line[0:2] == 'MT':
                spline = line.split("\t")
                spline[0] = "chrM"
                mitotemp.write("\t".join(spline))
            line = 'chr' + line
        temp.write(line)
    temp.close()
    mitotemp.close()
    os.system("cat mtmp tmp > "+hg18vcf+" ; rm mtmp ; rm tmp")

if __name__ == "__main__":
    main()
