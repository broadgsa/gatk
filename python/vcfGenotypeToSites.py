import sys
genotype_vcf = open(sys.argv[1])
sites_vcf = open(sys.argv[2],'w')

sample_name = "ALL_HET"
info = "."
format = "GT"
het = "0/1"
use_fields = range(7)

line_counter = 0
print("Reading genotype file...")
for line in genotype_vcf.readlines():
    line_counter += 1
    if ( line.startswith("#") and not line.startswith("#CHR") ):
        sites_vcf.write(line)
    elif ( line.startswith("#CHR") ):
        sites_vcf.write("##source=vcfGenotypeToSites\n")
        spline = line.strip().split("\t")
        newfields = list()
        for i in range(9):
            newfields.append(spline[i])
        newfields.append(sample_name)
        sites_vcf.write("\t".join(newfields)+"\n")
    else:
        spline = line.strip().split("\t")
        newfields = list()
        for i in use_fields:
            newfields.append(spline[i])
        newfields.append(info)
        newfields.append(format)
        newfields.append(het)
        sites_vcf.write("\t".join(newfields)+"\n")
    if ( line_counter % 100000 == 0 ):
        print("Converted: "+str(line_counter)+" lines")

genotype_vcf.close()
sites_vcf.close()
