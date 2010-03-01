#!/usr/bin/env python

import sys, FlatFileTable, os

if sys.argv < 3:
    print "Usage: AnnotateVCFwithMAF.py VCF_file MAF_file"
    sys.exit()
    
vcf_filename = sys.argv[1]
maf_filename = sys.argv[2]

maf_gen = FlatFileTable.record_generator(maf_filename, "\t")

headers=["gene","type","transcript","strand","genomechange","cDNAchange","codonchange","proteinchange"]

loci_and_info = []

for record in maf_gen:
    #print record
    #info_string = ",".join(["%s=%s" % (header, record[header]) for header in headers])
    info_string = ""
    for index,header in enumerate(headers):
        if record.has_key(header):
            if index > 0:
                info_string += ";"
            info_string += "%s=%s" % (header, record[header])

    locus = record["chr"]+":"+record["start"]
            
    #print locus, info_string
    loci_and_info.append((locus, info_string))
    
#vcf_gen = FlatFileTable.record_generator(vcf_file, "\t", 34)
vcf_file = open(vcf_filename)
vcf_out_file = open(os.path.splitext(os.path.basename(vcf_filename))[0]+".maf_annotated.vcf", "w")
vcf_format_line = vcf_file.readline()
vcf_out_file.write(vcf_format_line)
if vcf_format_line != "##fileformat=VCFv3.3\n":
    print ("VCF not v 3.3")
    sys.exit()

header = vcf_file.readline()
while header != "" and header.startswith("#"):
    if header.startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"):
        break
    vcf_out_file.write(header)
    header = vcf_file.readline()

header_fields = header
if not header_fields.startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"):
    print ("VCF header fields not in expected order")
    print header_fields
    sys.exit()

vcf_out_file.write("##source=AnnotateVCFwithMAF\n")
for header_field in headers:
    vcf_out_file.write("##INFO="+header_field+",1,String,"+header_field+"\n")
vcf_out_file.write(header_fields)

for vcf_line, locus_and_info in zip(vcf_file.readlines(), loci_and_info):
    vcf_line_fields = vcf_line.split("\t")
    vcf_locus = vcf_line_fields[0]+":"+vcf_line_fields[1]
    #print record
    maf_locus, maf_info = locus_and_info
    if maf_locus != vcf_locus:
        print "ERROR: VCF and MAF loci did not match"
        sys.exit()
    
    vcf_line_fields[7] = vcf_line_fields[7]+";"+maf_info
    new_vcf_line = "\t".join(vcf_line_fields)
    vcf_out_file.write(new_vcf_line)
