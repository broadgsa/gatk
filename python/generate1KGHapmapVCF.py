#!/usr/bin/env python

import os

hapmap_dir = os.getcwd()+"/" ##CHANGE ME

def convert(line):
    line = line.replace("chr","",1)
    if ( line.startswith("M") ):
        line = line.replace("M","MT",1)
    return line

for file in os.listdir(hapmap_dir):
    if ( file.endswith('vcf') ):
        chrM_lines = list()
        print("converting: "+file)
        in_vcf = open(hapmap_dir+file)
        out_vcf_filename = file.replace("hg18","b36")
        out_vcf = open(out_vcf_filename,'w')
        for line in in_vcf.readlines():
            if ( line.startswith("#") ):
                out_vcf.write(line)
            else:
                if ( line.startswith("chrM") ):
                    chrM_lines.append(line)
                else:
                    out_vcf.write(convert(line))
        for line in chrM_lines:
            out_vcf.write(convert(line))
