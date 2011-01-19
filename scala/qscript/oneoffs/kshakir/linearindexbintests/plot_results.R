#!/bin/env Rscript

require(lattice)
require(sqldf)

args <- commandArgs(TRUE)
verbose = TRUE

input = args[1]
memory = args[2]

mfpb_data <- read.table(input, head=T)
mfpb_data <- mfpb_data[order(mfpb_data$set, mfpb_data$max_features_per_bin, mfpb_data$memory_limit_gb, mfpb_data$run_number) , ]
mfpb_data <- sqldf("select \"set\", max_features_per_bin, memory_limit_gb, avg(cpu_s) as cpu_s, avg(max_memory_mb) as max_memory_mb from mfpb_data where job_success = 'done' group by \"set\", max_features_per_bin, memory_limit_gb")

outfile = paste("max_features_per_bin_Xmx", memory, "g.pdf", sep="")
pdf(outfile, height=7, width=14)

par(cex=1.3)

xyplot(max_memory_mb + cpu_s ~ log10(max_features_per_bin), groups = set, data = subset(mfpb_data, memory_limit_gb == memory), type="b", scales=list(relation="free"), auto.key=T)

dev.off()
