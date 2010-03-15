#!/usr/bin/env python

import sys
import os

bamfile_base = "/seq/picard_aggregation/"
fingerprint_base = "/seq/references/reference_genotypes/non-hapmap/"
hg18_reference = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"
hg18_dbsnp = "/humgen/gsa-hpprojects/GATK/data/dbsnp_130_hg18.rod"
b36_dbsnp = "/humgen/gsa-hpprojects/GATK/data/dbsnp_130_b36.rod"
b36_reference = "/broad/1KG/reference/human_b36_both.fasta"
hg18_intervals = "/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.targets.interval_list"
#hg18_intervals = "/humgen/gsa-hpprojects/FHS/indexed/interval_lists/fhs_jhs_pilot.targets.interval_list"
b36_intervals = ""

min_base_q = "10"
min_map_q = "10"
max_reads = "1000000"
min_conf = "50"
variant_expression = "QUAL <= 50.0 || AB > 0.75 || QD < 5.0 || HRun > 3"
spreadsheetPath = sys.argv[3]
projectName = sys.argv[2]
groupName = sys.argv[1]
reference = sys.argv[4]
filter_name = projectName+"_Initial_Filter"
if ( reference != "hg18" and reference != "b36" ):
    raise ValueError("Illegal reference type")
elif ( reference == "hg18" ):
    reference = hg18_reference
    dbsnp = hg18_dbsnp
    intervals = hg18_intervals
    fpref = "Homo_sapiens_assembly18"
else:
    reference = b36_reference
    dbsnp = b36_dbsnp
    intervals = b36_intervals
    fpref = "human_b36"

outputFile = projectName+"_bam_files.txt"
OUTPUT_HEADER = ["sample_id","recalibrated_bam_file","individual_id","fingerprint_file","reference_file","dbsnp_file","interval_list","max_reads_at_locus","min_confidence","min_mapping_quality","min_base_quality","variant_filter_expression","variant_filter_name"]

if ( spreadsheetPath.find("/") > -1 ):
    newSpreadsheet = spreadsheetPath.rsplit("/",1)[1].rsplit(".",1)[0]+"_proper_format.tsv"
else:
    newSpreadsheet = spreadsheetPath.rsplit(".",1)[0]+"_proper_format.tsv"

# convert to proper format
os.system("sed 's/\\r/\\n/g' "+spreadsheetPath+" > "+newSpreadsheet)

project_info = open(newSpreadsheet)

header = project_info.readline().strip().split("\t")

project_index = header.index("Project")
sample_index = header.index("Sample")
status_index = header.index("Sample Status")

def versionCompare(version1,version2):
    return -int(version1.split("v")[1])+int(version2.split("v")[1])

def getNewestVersion(baseDir):
    versions = os.listdir(baseDir)
    versions.sort(versionCompare)
    for version in versions:
        if ( "finished.txt" in os.listdir(baseDir+version+"/") ):
            return version

outputFile = open(outputFile,'w')
outputFile.write("\t".join(OUTPUT_HEADER)+"\n")

for line in project_info.readlines():
    if ( not line.startswith("\n") and not line.startswith(" ") and not line.startswith("\t") ):
        spline = line.strip().split("\t")
        versioningDirectory = bamfile_base+spline[project_index]+"/"+spline[sample_index]+"/"
        version = getNewestVersion(versioningDirectory)
        bamfile = versioningDirectory+version+"/"+spline[sample_index]+".bam"
        fingerprint_path = fingerprint_base+spline[project_index]+"/"+fpref+"/"
        if ( os.path.isdir(fingerprint_path) and spline[sample_index]+".fingerprint.geli" in os.listdir(fingerprint_path) ):
            fingerprint_file = fingerprint_path+spline[sample_index]+".fingerprint.geli"
        else:
            fingerprint_file = ""
        if ( spline[status_index] == "Complete" ):
            outputFile.write(projectName+"_"+spline[sample_index]+"\t"+bamfile+"\t"+groupName+"\t"+fingerprint_file+"\t"+reference+"\t"+dbsnp+"\t"+intervals+"\t"+max_reads+"\t"+min_conf+"\t"+min_map_q+"\t"+min_base_q+"\t"+variant_expression+"\t"+filter_name+"\n")
