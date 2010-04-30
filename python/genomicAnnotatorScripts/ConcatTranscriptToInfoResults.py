import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
from IndentedHelpFormatterWithNL import *

run_locally = True

# Init cmd-line args
description = """
This script creates and runs the command line that concatenates all 50 results of the GenerateTranscriptToInfo.py script into one big file that can be directly used with the GenomicAnnotator.
"""

parser = OptionParser( description=description, usage="usage: %prog [options] ", formatter=IndentedHelpFormatterWithNL())

parser.add_option("-r", "--refgene-directory", metavar="DIR", dest="refgene_dir", help="Specifies the directory that contains refGene-converted.txt", default="/humgen/gsa-hpprojects/GATK/data/Annotations/refseq/raw/")

(options, args) = parser.parse_args()

def error(msg):
    print("ERROR: %s.        (Rerun with -h to print help info) \n" % msg)
    #parser.print_help()
    sys.exit(-1)



contig_chars = ["M"] + range(1,23) + ["X", "Y"]

contigs = []
contigs += [ "chr" + str(x) for x in contig_chars ] 
contigs += [ "chr" + str(x) + "_random" for x in set( contig_chars ).difference(set(['M',12,14,20,'X','Y']))  ]    # There's no _random chromosomes for chrM,12,14,20,Y

#print(contigs)



# Update the refGene-big-table-header.txt header file using the header from one of the single-contig files.
command = "head -n 1 " + (options.refgene_dir + "/refGene-big-table-ucsc-%s.txt " % contigs[0]) + " > " + options.refgene_dir + "/refGene-big-table-header.txt"
print(command)
os.system(command)


# Concatenate
header_start = open(options.refgene_dir+"/refGene-big-table-header.txt").read().split("\t")[0]
command = "cat "
for contig in contigs:
    command += options.refgene_dir+"/refGene-big-table-ucsc-%s.txt " % contig
    
command += " | grep -v " + header_start 
command += " | cat " + options.refgene_dir+"/refGene-big-table-header.txt - "

command += " > " + options.refgene_dir+"/refGene-big-table-ucsc.txt"        
print(command)
os.system(command)

