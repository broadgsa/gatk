
#!/usr/bin/env python

import sys

dinuc_root = sys.argv[1]
fout = file(dinuc_root+".log_reg_params", "w")

fout.write("dinuc\t")
for p in range(5):
    for q in range(5):
        fout.write("%d,%d\t" % (p,q))
fout.write("\n")

for dinuc in ("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"):
    dinin = open(dinuc_root+"."+dinuc+".parameters")
    #dinin.readline()
    params = []
    for line in dinin:
        line.rstrip("\n")
        params.extend(map(float, line.split()))
    fout.write(dinuc+"\t")
    fout.write("\t".join(map(str, params)))
    fout.write("\n")
        







