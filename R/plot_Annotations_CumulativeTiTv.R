#!/broad/tools/apps/R-2.6.0/bin/Rscript

args <- commandArgs(TRUE)
verbose = TRUE

input = args[1]
annotationName = args[2]


c <- read.table(input, header=T)

#
# Plot residual error as a function of the covariate
#

gt = c[c$GT==1 & c$numVariants>1000,]
lt = c[c$GT==0 & c$numVariants>1000,]

outfile = paste(input, ".cumulativeTiTv_", annotationName, "_GTfilter.pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
plot(gt$value,gt$cumulativeTiTv,xlab=annotationName,ylab="Ti/Tv ratio",main=paste("Filter out SNPs with",annotationName,"> x",sep=" "),pch=20);
dev.off()

outfile = paste(input, ".cumulativeTiTv_", annotationName, "_LTfilter.pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
plot(lt$value,lt$cumulativeTiTv,xlab=annotationName,ylab="Ti/Tv ratio",main=paste("Filter out SNPs with",annotationName,"< x",sep=" "),pch=20);
dev.off()