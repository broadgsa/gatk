import sys
import os
import re
import traceback
import shlex, subprocess
from optparse import OptionParser, OptionGroup
from IndentedHelpFormatterWithNL import *
import time

# Init cmd-line args
description = """
This script submits LSF jobs that run the GATK TranscriptToInfo Walker on each individual chromosome. This reduces the overall runtime to a manageable ammount (eg. < 1 day).

NOTE: This script must be run in the top level dir of your GATK checkout area.
"""

parser = OptionParser( description=description, usage="usage: %prog [options] ", formatter=IndentedHelpFormatterWithNL())


parser.add_option("-i", "--transcript-table", metavar="PATH", dest="transcript_table", help="Path of the file that contains the transcript data in AnnotatorROD format (eg. /humgen/gsa-hpprojects/GATK/data/Annotations/refseq/raw/refGene-converted.txt)")
parser.add_option("-f", "--output-filename-prefix", metavar="PREFIX", dest="prefix", help="Output filename prefix (eg. refGene)")
parser.add_option("-p", "--print", dest="justprint", action="store_true", default=False, help="Only print the commands to standard out, don't actually execute them yet.")
parser.add_option("-e", "--execute", dest="execute", action="store_true", default=False, help="Executes the commands. This flag acts as a confirmation that you want to proceed with launching the processes.")
parser.add_option("-l", "--locally", dest="run_locally", action="store_true", default=False, help="Don't submit the commands to LSF. Run them sequentially on the current machine.")
parser.add_option("-R", "--reference", metavar="PATH", dest="reference", help="Specifies the path of the reference file to use.", default="/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
parser.add_option("-n", "--gene-name-columns", dest="gene_name_columns", metavar="GENE_NAMES", help="Comma-separated list of column names that contain gene names. This arg is passed through to the GenomicAnnotator. The GenomicAnnotator docs have more details on this.")
parser.add_option("-q", "--queue", dest="queue", metavar="QUEUE", help="Specifies the LSF queue to use.", default="solexa")
parser.add_option("-s", "--num-parallel-processes", dest="num_parallel_processes", metavar="SLOTS", help="How many processes to launch simultaneously. This is only used when the -l option is set.", default="1")

(options, args) = parser.parse_args()

def error(msg):
    print("ERROR: %s.        (Rerun with -h to print help info) \n" % msg)
    parser.print_help()
    sys.exit(-1)

run = options.execute
justprint = options.justprint
run_locally = options.run_locally

if not run and not justprint:
    error("Must run with either -p or -e")    

transcript_table = options.transcript_table
if not transcript_table or not os.access(transcript_table, os.R_OK):
    error("Must specify a valid transcript table file path using -i")    

gene_name_columns = options.gene_name_columns
if not gene_name_columns:
    error("Must specify gene name columns using -n")

output_file_prefix = options.prefix
if not output_file_prefix:
    error("Must specify the output file prefix using -f")

reference = options.reference
if not os.access(reference, os.R_OK):
    error("Couldn't access reference file: "+ reference)

queue = options.queue
num_parallel_processes = int(options.num_parallel_processes)

transcript_dir = os.path.dirname(transcript_table)
logs_dir = os.path.join(transcript_dir,"logs")

contig_chars = ["M"] + range(1,23) + ["X", "Y"]

contigs = []
contigs += [ "chr" + str(x) for x in contig_chars ] 
contigs += [ "chr" + str(x) + "_random" for x in set( contig_chars ).difference(set(['M',12,14,20,'X','Y']))  ]    # There are no "_random" chromosomes for chrM,12,14,20,Y




if run:
    print("Deleting any previous logs...")
    os.system("rm " + os.path.join(logs_dir,"bsub_*_log.txt"))
    os.system("mkdir " + logs_dir)


running_processes = []
def execute(command, stdout_filename=None):

    # Wait until a slot becomes open
    while len( running_processes ) + 1 >= num_parallel_processes:
        # Check if any have ended
        for process in running_processes:
            if process.poll() != None:
                print("Process [pid=" + str(process.pid) + "] finished with exit status: " + str(process.returncode))
                running_processes.remove(process)
                break
        else:
            time.sleep(3) # Sleep for 3 seconds before checking again

    # A slot has opened up - start another process
    stdout = None
    if stdout_filename:
        stdout = open(stdout_filename, "w+")
        

    p = subprocess.Popen(shlex.split(command), stdout=stdout, stderr=subprocess.STDOUT)
    running_processes.append(p)

for contig in contigs:        
    if contig.count("random") or contig.lower().count("chrm"):
        MEMORY_USAGE = 10  # Gigabytes
        EXCLUSIVE = ""
    else:
        if run_locally:
            MEMORY_USAGE = 64
        else:
            MEMORY_USAGE = 32
        EXCLUSIVE = ""
            
    command = "java -Xmx"+str(MEMORY_USAGE)+"g -jar dist/GenomeAnalysisTK.jar -T TranscriptToInfo -l info -nt " + str(num_parallel_processes) + " -R " + reference + " -B transcripts,AnnotatorInputTable,"+transcript_table+" -n "+gene_name_columns+" -o "+ os.path.join(transcript_dir,output_file_prefix) +"-big-table-ucsc-%s.txt -L %s:1+ " % (contig, contig)
    if not run_locally:
        command = "bsub "+EXCLUSIVE+" -q " + queue + " -R \"rusage[mem="+str(MEMORY_USAGE)+"]\" -o " + os.path.join(logs_dir,contig+"_log.txt") + " "  + command
    
    if run:
        print("Executing: " + command)
        if run_locally:
            #execute(command, os.path.join(logs_dir,contig+"_log.txt"))
            os.system(command + " >& " + os.path.join(logs_dir,contig+"_log.txt"))
        else:
            os.system(command)
    else:
        print(command)
