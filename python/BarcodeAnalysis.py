#!/usr/bin/python

import os,sys
from farm_commands import cmd

def get_line_number_from_file(filename, line_num):
    for index, line in enumerate(open(filename)):
        if index == line_num:
            return line

solexa_dir = "."

barcodes = []
barcode_file = open("../barcode_sequences.txt")
barcode_file.readline()
for line in barcode_file:
    fields = line.split()
    barcodes.append(fields[1])

farm = "long"
lane = 6
alias = "lane"+str(lane)
tmp_dir = "/humgen/gsa-scr1/andrewk/tmp/"

cmd("/seq/software/picard/current/bin/ExtractIlluminaBarcodes.jar B="+solexa_dir+" L="+str(lane)+" BARCODE_POSITION=1 METRICS_FILE=metrics.out BARCODE=" + " BARCODE=".join(barcodes))

cmd( "/seq/software/picard/current/bin/BustardToSam.jar B="+solexa_dir+" L="+str(lane)+" O="+alias+".bam SAMPLE_ALIAS="+alias+" RUN_BARCODE="+alias+" TMP_DIR="+tmp_dir)#, farm_queue=farm) #PAIRED_RUN=True 

for barcode in barcodes:
    id1 = cmd( "/seq/software/picard/current/bin/BustardToSam.jar B="+solexa_dir+" L="+str(lane)+" BARCODE_POSITION=1  BARCODE="+barcode+" O="+barcode+".bam RUN_BARCODE="+barcode+" SAMPLE_ALIAS=lane5 TMP_DIR="+tmp_dir, farm_queue=farm) #PAIRED_RUN=True 

    id2 = cmd( "/home/radon01/andrewk/dev/Sting/trunk/python/AlignBam.py "+barcode+".bam", farm_queue=farm, waitID = id1)
    id3 = cmd("/seq/software/picard/1.21/bin/MergeSamFiles.jar I="+barcode+".aligned.bam O="+barcode+".aligned.merged.bam VALIDATION_STRINGENCY=SILENT", farm_queue=farm, waitID = id2)
    id4 = cmd("java -Xmx4096m -jar /seq/software/picard/1.21/bin/MarkDuplicates.jar VALIDATION_STRINGENCY=SILENT I= "+barcode+".aligned.merged.bam  O="+barcode+".aligned.duplicates_marked.bam  M="+barcode+".aligned.bam.duplicates_marked.bam.duplicate_metrics TMP_DIR="+tmp_dir, farm_queue=farm), waitID=id3)
    cmd("java -jar /seq/software/picard/current/bin/CalculateHsMetrics.jar BI=/seq/references/HybSelOligos/thousand_genomes_alpha_redesign/thousand_genomes_alpha_redesign.baits.interval_list TI=/seq/references/HybSelOligos/thousand_genomes_alpha_redesign/thousand_genomes_alpha_redesign.targets.interval_list I="+barcode+".aligned.duplicates_marked.bam M="+barcode+".hs_metrics VALIDATION_STRINGENCY=SILENT", farm_queue=farm, waitID = id4)


fout = open("hybrid_selection.metrics", "w")
print >> fout, "BARCODE\t"+get_line_number_from_file(barcodes[0]+".hs_metrics", 6).rstrip()
for barcode in barcodes:
    print >> fout, barcode+"\t"+get_line_number_from_file(barcode+".hs_metrics", 7).rstrip()
