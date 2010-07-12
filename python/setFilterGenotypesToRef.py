import sys
print("Fixing "+sys.argv[1]+" to "+sys.argv[2])
bad_vcf = open(sys.argv[1])
out_vcf = open(sys.argv[2],'w')

for line in bad_vcf.readlines():
    if ( line.startswith("#") ):
        out_vcf.write(line)
    else:
        spline = line.strip().split("\t")
        newspline = list()
        for field in spline:
            if ( field.find("pGeno") > -1 ):
                field = "0/0:"+field.split(":",1)[1]
            newspline.append(field)
        out_vcf.write("\t".join(newspline)+"\n")
