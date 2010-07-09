#!/bin/env Rscript

args <- commandArgs(TRUE)
verbose = TRUE

input = args[1]

data = read.table(input,sep=",",head=T)
numCurves = (length(data) - 1)/3
maxSpec = max(data[,(1:numCurves)*3])

outfile = paste(input, ".variantROCCurve.pdf", sep="")
pdf(outfile, height=7, width=7)

par(cex=1.3)
plot(data$specificity1,data$sensitivity1, type="n", xlim=c(0,maxSpec),ylim=c(0,1),xlab="1 - Specificity",ylab="Sensitivity")
for(iii in 1:numCurves) {
	points(data[,iii*3],data[,(iii-1)*3+2],lwd=3,type="l",col=iii)
}
legend("bottomright", names(data)[(0:(numCurves-1))*3+1], col=1:numCurves,lwd=3)
dev.off()
