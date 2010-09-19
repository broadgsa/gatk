#!/usr/bin/python
import sys
input_file = sys.argv[1]
file_index = 1

headerLines = list()
intervals = open(input_file)
prevContig = None
outFile = None

def parseContig(line):
    if( line.find("-") > -1 ): ## format is chr:start-stop
        return line.split(":")[0]
    else:
        return line.split("\t")[0]

for line in open(input_file).readlines():
    if ( line.startswith("@") ):
        headerLines.append(line)
    else:
        thisContig = parseContig(line)
        if ( thisContig != prevContig ):
            file_index += 1
            try:
                newOutFile = open(sys.argv[file_index],'w')
                if ( outFile != None):
                    outFile.close()
                outFile = newOutFile
                for headerline in headerLines:
                    outFile.write(headerline)
            except IndexError:
                print("Error: fewer output files than contigs. Writing remainder to final file.")
            prevContig = thisContig
        outFile.write(line)
