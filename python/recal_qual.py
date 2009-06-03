#!/usr/bin/env python

# Hg18 ref and rods
reference = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"
dbsnp = "/humgen/gsa-scr1/GATK_Data/dbsnp.assembled_and_random_contigs_only.rod.out"

## 1KG ref and rods
reference = "/broad/1KG/reference/human_b36_both.fasta"
dbsnp = "/humgen/gsa-scr1/GATK_Data/dbsnp.1kg.rod.out"

# Flags - many of these are not necessary now but were used for debugging
training_range_str = "" #" -L chr1"
gatk_flags = ""
downsample_fraction = "" #" --DOWNSAMPLE_FRACTION 0.01"
read_group = "" #" --READ_GROUP SRR005814"
min_mapping_quality = " --MIN_MAPPING_QUALITY 1"
recal_limit = 0

javaexe = "java -ea"

justPrintCommands = False

import os, sys
sys.path.append(".")
if len(sys.argv) > 1:
    #from recal_args import *
    #print "Using args imported from recal_args.py"
#else:
    #from default_broad_recal_args import *
    #print "Using default Broad recalibration arguments"
    file = sys.argv[1]
    fileroot = file
    if file.endswith(".bam"):
        fileroot = fileroot.replace(".bam", "")
    fileroot = os.path.basename(fileroot)
    print "fileroot: "+fileroot
    #training_range_str = " -L chr1:1-20000000"

plot = False
recal_only = False
for arg in sys.argv:
    if arg == "-plot":
        plot = True
    if arg == "-recal_only":
        recal_only = True

def cmd(cmdstr):
    print ">>> Running: "+cmdstr
    if (not justPrintCommands):
        ret_val = str(os.system(cmdstr))
        print "<<< Returned: "+str(ret_val)

        if int(ret_val) != 0:
            sys.exit("Last command failed: "+cmdstr)

# Assure that BAM index is available
def assert_bam_index(bam_filename):
    if not os.path.exists(bam_filename+".bai"):
        cmd("samtools index "+bam_filename)

def ensure_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

# Process arguments
if recal_limit == 0:
    recal_limit_str = ""
else:
    recal_limit_str = "_"+str(recal_limit)

# This was the central directory where I put the files, I've changed it to "." 
# for simplicity now even though it probably needs to be read group aware
#covar_path = "/humgen/gsa-scr1/andrewk/recalibration_work"
covar_path = "."
central_fileroot= covar_path+"/"+fileroot

bam = central_fileroot+".bam"
if not os.path.lexists(bam):
    os.symlink(file, bam)

correction_str = ".corrected"+recal_limit_str
corrected_fileroot = central_fileroot+correction_str
corrected_bam = corrected_fileroot+".bam"

assert_bam_index(bam)
cmd(javaexe+" -jar ~andrewk/dev/Sting/trunk/dist/GenomeAnalysisTK.jar -T CountCovariates  -R "+reference+" --DBSNP "+dbsnp+" -mqs 40 -I "+bam+" --OUTPUT_FILEROOT "+central_fileroot+" -l INFO --CREATE_TRAINING_DATA"+training_range_str+downsample_fraction+gatk_flags+read_group+min_mapping_quality)

# Make plots
if plot:
    cmd("~andrewk/covariates/plot_q_emp_stated_hst.R "+central_fileroot+".empirical_v_reported_quality.csv")
    cmd("~andrewk/covariates/plot_qual_diff_v_cycle_dinuc.R "+central_fileroot)
    
cmd("python ~andrewk/dev/Sting/trunk/python/LogRegression.py "+central_fileroot)

if recal_only: # Stop if we are only making recalibration tables
    sys.exit()

cmd(javaexe+" -jar ~/dev/Sting/trunk/dist/GenomeAnalysisTK.jar -T LogisticRecalibration -I "+bam+" -R "+reference+" --logisticparamsfile "+central_fileroot+".recalibration_table --outputbamfile "+corrected_bam+" -l DEBUG -M "+str(recal_limit)+gatk_flags)
cmd("samtools index "+corrected_bam)    

cmd(javaexe+" -jar ~andrewk/dev/Sting/trunk/dist/GenomeAnalysisTK.jar -T CountCovariates  -R "+reference+" --DBSNP "+dbsnp+" -mqs 40 -I "+corrected_bam+" --OUTPUT_FILEROOT "+corrected_fileroot+" -l INFO"+gatk_flags+read_group+min_mapping_quality)

# Make plots
if plot:
    cmd("~andrewk/covariates/plot_q_emp_stated_hst.R "+corrected_fileroot+".empirical_v_reported_quality.csv")
    cmd("~andrewk/covariates/plot_qual_diff_v_cycle_dinuc.R "+corrected_fileroot)

# Merge corrected and uncorrected files
output_dir = "comparison_png"
ensure_dir(output_dir)
corrected_pngs = [f for f in os.listdir(".") if f.endswith(".png") and "corrected" in f]
print "\n".join(corrected_pngs)

for corrected_png in corrected_pngs:
    raw_png = corrected_png.replace(correction_str, "")
    both_png = corrected_png.replace(correction_str, "raw_and_corrected")
    cmd(" ".join(["montage", raw_png, corrected_png, "-geometry +0+0", output_dir+"/"+both_png]))



