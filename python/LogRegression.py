#!/usr/bin/env python

import subprocess, time, sys

dinucs = ("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")

#covar_work_dir = "/humgen/gsa-scr1/andrewk/recalibration_work"
#covar_table_dir = "/humgen/gsa-scr1/andrewk/recalibration_tables"
fileroot = sys.argv[1]
covar_counts_root = fileroot+".covariate_counts"
parameter_root = fileroot+".log_reg"
recal_table = fileroot+".recalibration_table"

# Run all R process to do regression and wait for completion
kids = []
for dinuc in dinucs:
    kids.append(subprocess.Popen(["/home/radon01/andrewk/covariates/logistic_regression.R", covar_counts_root, parameter_root, dinuc]))
    #kids.append(subprocess.Popen(["sleep", "8"]))

while any(kid.poll() is None for kid in kids):
    time.sleep(0.25)

fout = file(recal_table, "w")

fout.write("dinuc_v.1\t")
for p in range(5):
    for q in range(5):
        fout.write("%d,%d\t" % (p,q))
fout.write("\n")

for dinuc in dinucs:
    dinin = open(parameter_root+"."+dinuc+".parameters")
    #dinin.readline()
    params = []
    for line in dinin:
        line.rstrip("\n")
        params.extend(map(float, line.split()))
    fout.write(dinuc+"\t")
    fout.write("\t".join(map(str, params)))
    fout.write("\n")
        







