#!/bin/env Rscript

args <- commandArgs(TRUE)
verbose = TRUE

d = read.table(args[1],head=T)
outfile = args[2]
title = args[3]

# -----------------------------------------------------------------------------------------------
# plot timing
# -----------------------------------------------------------------------------------------------
pdf(outfile, height=5, width=8)
boxplot(d$walltime ~ d$operation, ylab = "Elapsed wall time in seconds [Log10 Scale]", log="y", main=title, cex.axis=0.75)
dev.off()
