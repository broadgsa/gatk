#!/usr/bin/env python
#import farm_commands
import os.path
import sys
import re
import time
import math

def grepFile(string, file):
    return grep(string, file.readlines())

def grep(string, list):
    expr = re.compile(string)
    return [elem for elem in list if expr.match(elem)]

def dictFromCombinedErrorCoverageFile(f):
    dictList = []
    for line in f:
        ln = line.split()
        key = ln[0]
        item = ln
        dictList.append([key,item])
    return dict(dictList)

def qualFromLod(L):
    X = math.exp(-L)
    try:
        return math.floor(-10*math.log10(X/1+X))
    except OverflowError:
        return 1000

callFileStr = sys.argv[1]
callFile = open(callFileStr)
pipelineSamplesFile = open(sys.argv[2])
directory = os.getcwd()
syzygyPathSamples = "/humgen/gsa-hphome1/flannick/pfizer/pspipeline/output/samples/"
sortByRefPath = "/humgen/gsa-hphome1/chartl/sting/perl/sortByRef.pl"

poolNames = []
poolInternalIDs = []
poolCombinedErrorCoverageCalls = []
proj = "" # required for an output file
for line in pipelineSamplesFile:
    ln = line.strip('\n')
    ln = ln.split(";")
    #ln = line.split("\t") # -depends on format
    poolNames.append(ln.pop(2))
    piid = ln.pop(0)
    poolInternalIDs.append(piid)
    proj = ln.pop(0)
    ceccPath = syzygyPathSamples+proj+"/"+piid+"/"+proj+"."+piid+".bam.combined.error.coverage.calls"
    print("reading: "+ceccPath)
    poolCombinedErrorCoverageCalls.append(dictFromCombinedErrorCoverageFile(open(ceccPath)))

pooledOutputFiles = []
for pool in poolNames:
    pooledOutputFiles.append(open(directory+"/"+pool+"_calls.vcf",'w'))

pooledCallsFile = open(directory+"/"+proj+"_combined_calls.vcf",'w')

header1 = "##format=PCFv1\n"
header2 = "##filedate="+time.strftime("%Y%m%d")+"\n"
header3 = "##source=expandedSummaryToVCF:"+callFileStr+"\n"
header4 = "##reference=Homo_sapiens_assembly18\n"
header5 = "##phasing=pooled\n"
header6 = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"+"\t"+"\t".join(poolNames)+"\n" #note: really cool method

# print header
pooledCallsFile.write(header1)
pooledCallsFile.write(header2)
pooledCallsFile.write(header3)
pooledCallsFile.write(header4)
pooledCallsFile.write(header5)
pooledCallsFile.write(header6)

# print rest of headers

for i in range(len(poolNames)):
    header6 = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    pooledOutputFiles[i].write(header1)
    pooledOutputFiles[i].write(header2)
    pooledOutputFiles[i].write(header3)
    pooledOutputFiles[i].write(header4)
    pooledOutputFiles[i].write(header5)
    pooledOutputFiles[i].write(header6)

for line in callFile:
# note: file is a .csv; so comma delimited with headers
# chr:position,alleles,gene,type,rs#,base_chg,annot,dist,pop_nr_freq,f_nr_case,f_nr_con,chisqstat_assoc,lrt_stat,f+,f-,SLOD,p_val,min_snp_dist,min_amp_dist,num_sig_pools,mean_worst_dir_nonref, ...
# ..., max_worst_dir_nonref,mean_diff_dir_nonref,max_diff_dir_nonref,mean_other/minor,max_other/minor,mean_frac_non_concord,max_frac_non_concord,mean_combined_lod,max_combined_lod,mean_cov_tot, ...
# ..., max_cov_tot,mean_cov_diff,max_cov_diff,mean_lod,max_lod,mean_lod_diff,max_lod_diff,mean_lod_diff_norm,max_lod_diff_norm,target_name,sig_pools

# call file assumed already to have been split by pool

# pileupFile is a .bam.coverage file from the pooled pipeline

# call file isn't sorted by chr:position 
# make the header strings

    if line.startswith("chr:pos"):
        continue
    elif not line.startswith("chr"):
        continue
    else:
        #print(line)
        entries = line.split(",")
        chrompos = entries[0]
        alleles = entries[1]
        variant = alleles[1]
        ref = alleles[0]
        dbSNP = entries[4]
        if dbSNP:
            pass
        else:
            dbSNP="."

        supportingPools = entries.pop(62).rstrip(']').lstrip('[')
        supportingPools = supportingPools.split(";")
        total_quality = 0
        total_slod = 0
        total_depth = 0
        quality_by_pool = []
        slod_by_pool = []
        depth_by_pool = []
        #sys.exit()
        for i in range(len(poolNames)):
            #grab line from the correct dict
            depth = 0;
            slod = 0;
            qual = 0;
            try:
                ceccLine = poolCombinedErrorCoverageCalls[i][chrompos]
                depth = ceccLine[18]
                qual=qualFromLod(float(ceccLine[21]))
                if ceccLine[22] == "NA":
                    slod = 0;
                else:
                    try:
                        slod = math.log10(float(ceccLine[23]))
                    except OverflowError:
                        slod = -1000;
            except KeyError:
                # do nothing
                pass
            #print this out to the file
            chromsplit = chrompos.split(":")
            outstr=chromsplit[0]+"\t"+chromsplit[1]+"\t"+dbSNP+"\t"+ref+"\t"+variant+"\t"+str(qual)+"\t0\t"+"DP="+str(depth)+";SB="+str(slod)+"\n"
            pooledOutputFiles[i].write(outstr)
            #now update data
            total_slod = total_slod + float(slod)
            total_depth = total_depth + int(depth)
            total_quality = total_quality + qual
            depth_by_pool.append(depth)
            slod_by_pool.append(slod)
            quality_by_pool.append(qual)
        #now for the pooled file
        chromsplit=chrompos.split(":")
        outstr = chromsplit[0]+"\t"+chromsplit[1]+"\t"+dbSNP+"\t"+ref+"\t"+variant+"\t"+str(total_quality)+"\t0\t"+"DP="+str(total_depth)+";SB="+str(total_slod)+";NP="+str(len(supportingPools))+"\tGT:GQ:DP:SB"
        pooledCallsFile.write(outstr)
        #propagate individual pool information
        for i in range(len(poolNames)):
            phase = "0/0"
            if grep(poolNames[i],supportingPools):
                phase = "0/1"
            else:
                phase = "0/0"
            pooledOut="\t"+phase+":"+str(quality_by_pool[i])+":"+str(depth_by_pool[i])+":"+str(slod_by_pool[i])
            pooledCallsFile.write(pooledOut)
        pooledCallsFile.write("\n")

## close all files ##

pooledCallsFile.close()

for i in range(len(poolNames)):
    pooledOutputFiles[i].close()

## sort the files -- system commands ##

for pool in poolNames:
    cmd1 = "cat "+directory+"/"+pool+"_calls.vcf | sort -n -k2,2 | perl "+sortByRefPath+" - /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta.fai > "+directory+"/tmp.vcf"
    cmd2 = "cp "+directory+"/tmp.vcf "+directory+"/"+pool+"_calls.vcf"
    os.system(cmd1)
    os.system(cmd2)

cmd = "cat "+directory+"/"+proj+"_combined_calls.vcf | sort -n -k2,2 | perl "+sortByRefPath+" - /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta.fai > "+directory+"/tmp.vcf"
os.system(cmd)
cmd = "cp "+directory+"/tmp.vcf "+directory+"/"+proj+"_combined_calls.vcf"
os.system(cmd)
cmd = "rm "+directory+"/tmp.vcf"
os.system(cmd)
