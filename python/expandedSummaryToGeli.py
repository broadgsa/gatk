#!/usr/bin/env python
import farm_commands
import os.path
import sys
import re


def grepFile(string, file):
    return grep(string, file.readlines())

def grep(string, list):
    expr = re.compile(string)
    return [elem for elem in list if expr.match(elem)]

callFile = open(sys.argv[1])
pileupFileName = sys.argv[2]

for line in callFile:
# note: file is a .csv; so comma delimited with headers
# chr:position,alleles,gene,type,rs#,base_chg,annot,dist,pop_nr_freq,f_nr_case,f_nr_con,chisqstat_assoc,lrt_stat,f+,f-,SLOD,p_val,min_snp_dist,min_amp_dist,num_sig_pools,mean_worst_dir_nonref, ...
# ..., max_worst_dir_nonref,mean_diff_dir_nonref,max_diff_dir_nonref,mean_other/minor,max_other/minor,mean_frac_non_concord,max_frac_non_concord,mean_combined_lod,max_combined_lod,mean_cov_tot, ...
# ..., max_cov_tot,mean_cov_diff,max_cov_diff,mean_lod,max_lod,mean_lod_diff,max_lod_diff,mean_lod_diff_norm,max_lod_diff_norm,target_name,sig_pools

# call file assumed already to have been split by pool

# pileupFile is a .bam.coverage file from the pooled pipeline

# call file isn't sorted by chr:position 
    if line.startswith("chr:pos"):
        continue
    elif not line.startswith("chr"):
        continue
    else:
        #print(line)
        entries = line.split(",")
        chrompos = entries.pop(0) # (now 'alleles' is 0-element)
        variant = entries.pop(0) # (now 'gene' is 0-element)
        #print(chrompos)
        #print(variant)
        #change format of "chr:pos" to "chr pos" for lookup
        #somehow this next line don't work
        #chrompos.replace(":"," ")
        g = chrompos.split(":");
        chrompos = g.pop(0)+" "+g.pop(0)
        #print(chrompos)
        pileupLine = grepFile(chrompos,open(pileupFileName))
        #print(pileupLine)
        # line is
        # chr pos ref num_A num_C num_G num_T
        if not pileupLine:
            continue
        else:
            pileupLine = pileupLine.pop(0)
            #print(pileupLine)
            pileupList = pileupLine.split(" ")
            ref = pileupList.pop(2) # now num_A is 2
            num_A = int(pileupList.pop(2))
            num_C = int(pileupList.pop(2))
            num_G = int(pileupList.pop(2))
            num_T = int(pileupList.pop(2))
            depth = num_A+num_C+num_G+num_T
            # output is
            # chr pos ref depth mapping call btr btnb AA AC AG AT CC CG CT GG GT TT
            outStr = chrompos+" "+ref+" "+str(depth)+" 5 "+variant+" -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12"
            print(outStr)
