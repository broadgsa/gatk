#!/usr/bin/env python
#import farm_commands
import subprocess
import os
import sys
import re
import time
import math

dbsnp = "."  # can be global. Let something else annotate dbsnp info
filter = "0" # don't declare any filtering in the VCF file

#print(sys.argv)

raw_calls_file = open(sys.argv[1])
output_vcf_file = open(sys.argv[2],'w')
pool_name = []

try:
    pool_name = sys.argv[3]
except IndexError:
    # parse the file name
    filepath = sys.argv[0].strip().split("/")
    pool_name = filepath[len(filepath)]


header = raw_calls_file.readline().strip()

fields = header.split()

#print(fields)

# parse it for important offsets

chrompos = fields.index("chr:offset")
ref_offset = fields.index("ref_base")
for_a_index = fields.index("A")
for_c_index = fields.index("C")
for_g_index = fields.index("G")
for_t_index = fields.index("T")
for_d_index = fields.index("D")
for_i_index = fields.index("I")
for_depth_index = fields.index("sum")
rev_a_index = fields.index("AR")
rev_c_index = fields.index("CR")
rev_g_index = fields.index("GR")
rev_t_index = fields.index("TR")
rev_d_index = fields.index("DR")
rev_i_index = fields.index("IR")
rev_depth_index = fields.index("sumr")
total_depth_index = fields.index("combined_sum")
lod_score_index = fields.index("combined_lod")
call_index = fields.index("allele2")
fish_flag = fields.index("flag")
fish_pval_index = fields.index("fisher-pval")

# now print the VCF header
head1 = "##source=Syzygy"
head2 = "##format=VCRv3.2"
fields = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",pool_name]
head3 = "#" + "\t".join(fields)

output_vcf_file.write(head1+"\n")
output_vcf_file.write(head2+"\n")
output_vcf_file.write(head3+"\n")

def getProportionNonref(list):
    total_bases = int(list[total_depth_index]) - int(list[rev_i_index].split(":")[1]) - int(list[rev_d_index].split(":")[1]) - int(list[for_i_index].split(":")[1]) - int(list[for_d_index].split(":")[1])
    ref_base = list[ref_offset]
    if ( ref_base == "A" ):
        ref_bases = int(list[for_a_index].split(":")[1]) + int(list[rev_a_index].split(":")[1])
    elif ( ref_base == "C"):
        ref_bases = int(list[for_c_index].split(":")[1]) + int(list[rev_c_index].split(":")[1])
    elif ( ref_base == "G"):
        ref_bases = int(list[for_g_index].split(":")[1]) + int(list[rev_g_index].split(":")[1])
    else:
        ref_bases = int(list[for_t_index].split(":")[1]) + int(list[rev_t_index].split(":")[1])

    return 1.0 - ( float(ref_bases+1) / float(total_bases+1) )

def generateVCFLine(chrom, pos, db, ref, alt, filt, qual, INFO):
    # make the info into a single string
    infoString = ""
    for keval in INFO:
        info = "=".join(keval)
        infoString = infoString+info
    format = "GT:GQ"
    genotype = "0/1:"+qual
    all_fields = [chrom, pos, db, ref, alt, qual, filt, infoString, format, genotype]
    return "\t".join(all_fields)

# instantiate a line buffer

two_lines_ago = [];
previous_line = [];
this_line = [];
next_line = [];
line_after_next = [];

# read through the file

for line in raw_calls_file.readlines():
    spline = line.strip().split() 

    # iterate through the lines

    two_lines_ago = previous_line
    previous_line = this_line
    this_line = next_line
    next_line = line_after_next
    line_after_next = spline

    # window has been updated

    if ( not two_lines_ago ):
        continue # window not filled yet
    else:
        if ( float(this_line[lod_score_index]) > 0 and this_line[call_index] != "D" and this_line[call_index] != "I"):
            # potential call
            chrom = this_line[chrompos].split(":")[0]
            pos = this_line[chrompos].split(":")[1]
            ref = this_line[ref_offset]
            alt = this_line[call_index]
            lod = this_line[lod_score_index]
            
            # standard vcf info made -- now add INFO
            # syzy depth
            syz_depth = this_line[total_depth_index]
            # syzy strand bias
            syz_sb = ""
            if ( this_line[fish_flag] == "NA" ):
                syz_sb = "-1"
            else:
                syz_sb = this_line[fish_pval_index]

            # do we want any other kind of annotations here ??

            # syzy neighborhood mismatch rate
            next_mmr = getProportionNonref(next_line) # from later in genome
            after_next_mmr = getProportionNonref(line_after_next) # from same
            prev_mmr = getProportionNonref(previous_line) # from earlier
            before_last_mmr = getProportionNonref(two_lines_ago) # same
            syzy_nmmr = str(next_mmr+after_next_mmr+prev_mmr+before_last_mmr)
            # turn these into key value pairs
            INFO = [["syzy_DP",syz_depth],[";syzy_SB",syz_sb],[";syzy_NMMR",syzy_nmmr]] # semicolons are a complete hack
            # get the vcf line
            vcfLine = generateVCFLine(chrom, pos, dbsnp, ref, alt, filter, lod, INFO)
            # print the sucker
            output_vcf_file.write(vcfLine+"\n")
