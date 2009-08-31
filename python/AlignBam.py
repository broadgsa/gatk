#!/usr/bin/env python

import os,sys
from farm_commands import cmd

def run_paired_BWA(bam):
    "Takes a paired end BAM file and outputs an aligned '.aligned.bam' file "
    lane = 0
    root = os.path.splitext(bam)[0]
    
    # BAM to BFQ
    cmd("/seq/software/picard/1.21/bin/BamToBfq.jar I="+bam+" LANE="+str(lane)+" FLOWCELL_BARCODE="+root+" PAIRED_RUN=True RUN_BARCODE="+root+" ANALYSIS_DIR=.")
    
    # BFQ to FQ
    root1 = root+"."+str(lane)+".0.1"
    root2 = root+"."+str(lane)+".0.2"
    bfq1 = root1+".bfq"
    bfq2 = root2+".bfq"
    fq1 = root1+".fq"
    fq2 = root2+".fq"

    cmd("/seq/dirseq/maq-0.7.1/maq bfq2fastq "+bfq1+" "+fq1)
    cmd("/seq/dirseq/maq-0.7.1/maq bfq2fastq "+bfq2+" "+fq2)

    # Run BWA
    sai1 = root1+".sai"
    sai2 = root2+".sai"
    cmd("/seq/dirseq/bwa/bwa-0.4.9/bwa aln /humgen/gsa-scr1/ebanks/bwaref/Homo_sapiens_assembly18.fasta "+fq1+" >"+sai1)
    cmd("/seq/dirseq/bwa/bwa-0.4.9/bwa aln /humgen/gsa-scr1/ebanks/bwaref/Homo_sapiens_assembly18.fasta "+fq2+" >"+sai2)
    
    # Pairing
    alignedsam = root+".aligned.sam"
    cmd("/seq/dirseq/bwa/bwa-0.4.9/bwa sampe /humgen/gsa-scr1/ebanks/bwaref/Homo_sapiens_assembly18.fasta "+sai1+" "+sai2+" "+fq1+" "+fq2+" > "+alignedsam)
    
    # SAM to BAM
    bam_ref_list = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.bam.ref_list"
    alignedbam = root+".aligned.bam"
    cmd("/seq/dirseq/samtools/current/samtools import "+bam_ref_list+" "+alignedsam+" "+alignedbam)
    
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: AlignBam.py BAM_FILE"
        sys.exit()
    bam = sys.argv[1]
    run_paired_BWA(bam)
    
        
