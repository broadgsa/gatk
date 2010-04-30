import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
from IndentedHelpFormatterWithNL import *

run_locally = True

# Init cmd-line args
description = """
This script runs a command that concatenates all 50 results of the GenerateTranscriptToInfo.py script into one big file that can be directly used by the GenomicAnnotator.
"""

parser = OptionParser( description=description, usage="usage: %prog [options] ", formatter=IndentedHelpFormatterWithNL())

parser.add_option("-r", "--refgene-directory", metavar="DIR", dest="refgene_dir", help="Specifies the directory that contains refGene-converted.txt", default="/humgen/gsa-hpprojects/GATK/data/Annotations/refseq/raw/")

parser.add_option("-u", "--ucsc", dest="ucsc", action="store_true", default=False, help="Generate the output file for use with the NCBI reference genome (this effects chromosome order and naming (eg. M chromosome is first and its called 'chrM' instead of 'MT')).")
parser.add_option("-n", "--ncbi", dest="ncbi", action="store_true", default=False, help="Generate the output file for use with the UCSC reference genome (this effects chromosome order and naming (eg. MT chromosome is last and its called 'MT' instead of 'chrM')).")

(options, args) = parser.parse_args()


def error(msg):
    print("ERROR: %s.        (Rerun with -h to print help info) \n" % msg)
    parser.print_help()
    sys.exit(-1)

ucsc = options.ucsc
ncbi = options.ncbi

if not ucsc and not ncbi: 
    error("Must run with either -u or -n")    


contig_chars = []
if ucsc:
    contig_chars = ["M"] + range(1,23) + ["X", "Y"]    
else:
    contig_chars = range(1,23) + ["X", "Y", "M"]


contigs = []
contigs += [ "chr" + str(x) for x in contig_chars ] 

if ucsc:  # NCBI doesn't have the _random contigs
    contigs += [ "chr" + str(x) + "_random" for x in set( contig_chars ).difference(set(['M','MT',12,14,20,'X','Y']))  ]    # There's no _random chromosomes for chrM,12,14,20,Y

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
if ncbi:
    command += "| perl -pe 's/^chrM(.*)$/MT\\1/i' | perl -pe 's/^chr([^p].*)$/\\1/i' "    # rename chrM to MT and remove the 'chr' from chromosome names 

command += " | cat " + options.refgene_dir+"/refGene-big-table-header.txt - "

if ucsc:
    command += " > " + options.refgene_dir+"/refGene-big-table-ucsc.txt"        
else:
    command += " > " + options.refgene_dir+"/refGene-big-table-ncbi.txt"        

print(command)
os.system(command)

