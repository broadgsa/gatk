#!/bin/env Rscript

args <- commandArgs(TRUE)
verbose = TRUE

input = args[1]
annotationName = args[2]

data = read.table(input,sep=",",head=T)

outfile = paste(input, ".ClusterReport.pdf", sep="")
pdf(outfile, height=7, width=8)

maxP = max(data$knownDist, data$novelDist)

plot(data$annotationValue, data$knownDist, ylim=c(0,maxP),type="b",col="orange",lwd=2,xlab=annotationName,ylab="fraction of SNPs")
points(data$annotationValue, data$novelDist, type="b",col="blue",lwd=2)
legend('topright', c('knowns','novels'),lwd=2,col=c("orange","blue"))
dev.off()
