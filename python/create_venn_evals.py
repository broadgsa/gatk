#!/usr/bin/env python

# Creates summary stats from a concordance file - SNPs called and number in dbSNP - and will output 

import sys, re, os, farm_commands

sting_dir = os.environ.get("STING_DIR")
if sting_dir == None:
    print "You need to set the environment variable STING_DIR to point to your Sting/GATK build"
    sys.exit()
    
if len(sys.argv) < 2:
    print "Usage: PROG CONCORDANCE_FILE"
    sys.exit()

class callset_data:
    def __init__(self, snps, in_dbSNP, titv):
        self.snps = snps
        self.in_dbSNP = in_dbSNP
        self.titv = titv
    def inc_snps(self):
        self.snps += 1
    def inc_dbsnps(self):
        self.in_dbSNP += 1
    def __str__(self):
        return "SNPs: %10d  in_dbSNP: %d  TiTv: %d" % (self.snps, self.in_dbSNP, self.titv)

sets = dict()
concordance_filename = os.path.abspath(sys.argv[1])
for line in file(concordance_filename):
    fields = line.split("\t")
    filter = 0
    if not re.search('HybridSelectionVariantFilter', line):
        match = re.search('Venn=(\S*)', line) ### TODO: Only works in 2-way venn with ";", remove ";" for N-way Venn - should be changed to work for both automatically
        if match:
            my_set = match.groups()[0]
            
            if sets.has_key(my_set):
                sets[my_set].inc_snps()
                if fields[2] != ".":
                    sets[my_set].inc_dbsnps()
            else:
                sets[my_set] = callset_data(1, 0, 0)

print "\n".join(["%40s %s" % (k,str(v)) for k,v in sets.items()])

print "Splitting concordance file\n"
splits_dir = concordance_filename+".splits"

for vennSet in sets.keys():
    print vennSet
    output_head = splits_dir+"_"+vennSet
    vcf_file = splits_dir+"_"+vennSet+".vcf"
    varianteval_file = splits_dir+"_"+vennSet+".varianteval"

    if not os.path.exists(varianteval_file):
        fout = open(vcf_file, "w")
        for line in file(concordance_filename):    
            if line.startswith("#"):
                fout.write(line)
            else:
                match = re.search('Venn='+vennSet, line)
                if match:
                    fout.write(line)
                    
        varianteval_out = os.path.basename(concordance_filename)+".vcf"
        farm_commands.cmd("java -Xmx4096m -jar "+sting_dir+" -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta "
                  +"-T VariantEval "
                  +"-B eval,VCF,"+vcf_file+" "
                  +"-o "+varianteval_file+" "
                  +"-D /humgen/gsa-scr1/GATK_Data/dbsnp_129_hg18.rod",
                  "gsa")









