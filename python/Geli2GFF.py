#!/usr/bin/env python

# Converts Geli files (genotype likelihood in picard) to GFF format for GATK

import sys, os, subprocess

Index2AffyChr = []
Index2AffyChr.append("chrM")
for i in range(1,23):
    Index2AffyChr.append("chr"+str(i))
Index2AffyChr.append("chrX")
Index2AffyChr.append("chrY")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "usage: Geli2GFF.py INPUT_FILE OUTPUT_DIRECTORY"
        sys.exit()
    else:
        infile = sys.argv[1]
        outdir = sys.argv[2]

    output_file = outdir+"/"+os.path.splitext(os.path.basename(infile))[0]+".gff"
    outfile = open(output_file, "w")

    cmd = "/seq/software/picard/current/bin/GeliToText.jar I="+infile
    pipe = subprocess.Popen(cmd, shell=True, bufsize=4000, stdout=subprocess.PIPE).stdout

    pipe.readline()
    for line in pipe:
        if line[0] not in ('@', '#'):
            fields = line.split("\t")
            try:
                contig = Index2AffyChr[int(fields[0])]
            except KeyError, ValueError:
                print "skipping "+fields[0]
                continue # Skip entries not in chr 0 through 24
                                      
            print >>outfile, contig+"\tHapmapGenotype\t"+fields[5]+"\t"+fields[1]+"\t"+str(int(fields[1])+1)+"\t.\t.\t.\t"


            






    
