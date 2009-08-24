#!/usr/bin/python

import os,sys
from farm_commands import cmd

def get_line_number_from_file(filename, line_num):
    for index, line in enumerate(open(filename)):
        if index == line_num:
            return line

if __name__ == "__main__":
    
    from optparse import OptionParser
    class OptionParser (OptionParser):
        def check_required (self, opts):
            for opt in opts:
                option = self.get_option(opt)
    
                # Assumes the option's 'default' is set to None!
                if getattr(self.values, option.dest) is None:
                    self.error("%s option not supplied" % option)
    
    parser = OptionParser()
    parser.add_option("-q", "--queue", help="Queue to submit jobs to (default: long)", dest="queue", default="long")
    parser.add_option("-l", "--lane", help="Lane number to process", dest="lane", default="")
    parser.add_option("-t", "--tmpdir", help="Temp directory to use", dest="tmpdir")
    parser.add_option("-b", "--barcodes", help="File with barcodes as \"Id[TABa]barcode_sequence\" with one header line", dest="barcodes")
    parser.add_option("-a", "--baits", help="Bait interval list (default: 1KG bait list)", dest="ti", default="/seq/references/HybSelOligos/thousand_genomes_alpha_redesign/thousand_genomes_alpha_redesign.baits.interval_list")
    parser.add_option("-r", "--targets", help="Target interval list (default: 1KG target list)", dest="ti", default="/seq/references/HybSelOligos/thousand_genomes_alpha_redesign/thousand_genomes_alpha_redesign.targets.interval_list")
    (options, args) = parser.parse_args()
    parser.check_required(("-l","-t","-b"))
        
    solexa_dir = "."
    
    barcodes = []
    barcode_file = open(options.barcodes) #"../barcode_sequences.txt")
    barcode_file.readline()
    for line in barcode_file:
        fields = line.split()
        barcodes.append(fields[1])
    
    farm = options.queue 
    lane = options.lane
    alias = "lane"+str(lane)
    tmp_dir = options.tmpdir #"/humgen/gsa-scr1/andrewk/tmp/"
    
    cmd("/seq/software/picard/current/bin/ExtractIlluminaBarcodes.jar B="+solexa_dir+" L="+str(lane)+" BARCODE_POSITION=1 METRICS_FILE=metrics.out BARCODE=" + " BARCODE=".join(barcodes))
    
    cmd( "/seq/software/picard/current/bin/BustardToSam.jar B="+solexa_dir+" L="+str(lane)+" O="+alias+".bam SAMPLE_ALIAS="+alias+" RUN_BARCODE="+alias+" TMP_DIR="+tmp_dir)#, farm_queue=farm) #PAIRED_RUN=True 
    
    for barcode in barcodes:
        id1 = cmd( "/seq/software/picard/current/bin/BustardToSam.jar B="+solexa_dir+" L="+str(lane)+" BARCODE_POSITION=1  BARCODE="+barcode+" O="+barcode+".bam RUN_BARCODE="+barcode+" SAMPLE_ALIAS=lane5 TMP_DIR="+tmp_dir, farm_queue=farm) #PAIRED_RUN=True 
    
        id2 = cmd( "/home/radon01/andrewk/dev/Sting/trunk/python/AlignBam.py "+barcode+".bam", farm_queue=farm, waitID = id1)
        id3 = cmd("/seq/software/picard/1.21/bin/MergeSamFiles.jar I="+barcode+".aligned.bam O="+barcode+".aligned.merged.bam VALIDATION_STRINGENCY=SILENT", farm_queue=farm, waitID = id2)
        id4 = cmd("java -Xmx4096m -jar /seq/software/picard/1.21/bin/MarkDuplicates.jar VALIDATION_STRINGENCY=SILENT I= "+barcode+".aligned.merged.bam  O="+barcode+".aligned.duplicates_marked.bam  M="+barcode+".aligned.bam.duplicates_marked.bam.duplicate_metrics TMP_DIR="+tmp_dir, farm_queue=farm, waitID=id3)
        cmd("java -jar /seq/software/picard/current/bin/CalculateHsMetrics.jar BI="+options.bait_list+" TI="+options.target_list+" I="+barcode+".aligned.duplicates_marked.bam M="+barcode+".hs_metrics VALIDATION_STRINGENCY=SILENT", farm_queue=farm, waitID = id4)
    
    
    fout = open("hybrid_selection.metrics", "w")
    print >> fout, "BARCODE\t"+get_line_number_from_file(barcodes[0]+".hs_metrics", 6).rstrip()
    for barcode in barcodes:
        print >> fout, barcode+"\t"+get_line_number_from_file(barcode+".hs_metrics", 7).rstrip()
    
