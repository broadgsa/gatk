import sys
import os
import re
import traceback
from optparse import OptionParser, OptionGroup
from IndentedHelpFormatterWithNL import *

run_locally = False

# Init cmd-line args
description = """
This script submits LSF jobs that run the GATK TranscriptToInfo Walker on each individual chromosome. This reduces the overall runtime to a managable ammount (eg. < 1 day).

NOTE: This script must be run in the top level dir of your GATK checkout area.
"""

parser = OptionParser( description=description, usage="usage: %prog [options] ", formatter=IndentedHelpFormatterWithNL())

parser.add_option("-d", "--refgene-directory", metavar="DIR", dest="refgene_dir", help="Specifies the directory that contains refGene-converted.txt", default="/humgen/gsa-hpprojects/GATK/data/Annotations/refseq/raw/")

parser.add_option("-p", "--print", dest="output", action="store_true", default=False, help="Only print the commands to standard out, don't actually execute them yet.")
parser.add_option("-e", "--execute", dest="execute", action="store_true", default=False, help="Executes the commands. This flag acts as a confirmation that you want to proceed with launching the processes.")
parser.add_option("-l", "--locally", dest="run_locally", action="store_true", default=False, help="Don't submit the commands to LSF. Run them sequentially on the current machine.")

(options, args) = parser.parse_args()

def error(msg):
    print("ERROR: %s.        (Rerun with -h to print help info) \n" % msg)
    parser.print_help()
    sys.exit(-1)

run = options.execute
output = options.output
run_locally = options.run_locally

if not run and not output: 
    error("Must run with either -p or -e")    




contig_chars = ["M"] + range(1,23) + ["X", "Y"]

contigs = []
contigs += [ "chr" + str(x) for x in contig_chars ] 
contigs += [ "chr" + str(x) + "_random" for x in set( contig_chars ).difference(set(['M',12,14,20,'X','Y']))  ]    # There are no "_random" chromosomes for chrM,12,14,20,Y

#print(contigs)


if run:
    print("Deleting any previous logs...")
    os.system("rm " + options.refgene_dir+"/logs/bsub_*_log.txt")
for contig in contigs:
        
    if contig.count("random") or contig.lower().count("chrm"):
        MEMORY_USAGE = 10  #Gigabytes
        EXCLUSIVE = ""
    else:
        if run_locally:
            MEMORY_USAGE = 32
        else:
            MEMORY_USAGE = 15
        EXCLUSIVE = ""
            
    command = "java -Xmx"+str(MEMORY_USAGE)+"g -jar dist/GenomeAnalysisTK.jar -T TranscriptToInfo -l info -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -B refgene,AnnotatorInputTable,"+options.refgene_dir+"/refGene-converted.txt -o "+options.refgene_dir+"/refGene-big-table-ucsc-%s.txt -L %s:1+ " % (contig, contig)
        #print(command)
    if not run_locally:
        command = "bsub "+EXCLUSIVE+" -q solexa -R \"rusage[mem="+str(MEMORY_USAGE)+"]\" -o "+options.refgene_dir+"/logs/bsub_"+contig+"_log.txt "+command
    
    
    if run:
        print("Executing: " + command)
        os.system(command)
    else:
        print(command)
