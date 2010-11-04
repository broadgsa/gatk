import os.path
import sys
from optparse import OptionParser


def main():
    global OPTIONS
    usage = "usage: %prog [options] cmap input.soapsnp"
    parser = OptionParser(usage=usage)
    
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("Requires exactly 2 arguments")

    cmap = dict([reversed(line.split()) for line in open(args[0])])
    #print cmap
    for line in open(args[1]):
        parts = line.split()
        mapped = cmap[parts[0]]
        print "\t".join([mapped] + parts[1:])

if __name__ == "__main__":
    main()