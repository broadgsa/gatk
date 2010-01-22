#!/broad/tools/apps/R-2.6.0/bin/Rscript

args <- commandArgs(TRUE)
verbose = TRUE

input = args[1]
outputDir = args[2]
annotationName = args[3]
minBinCutoff = as.numeric(args[4])


c <- read.table(input, header=T)

#
# Plot TiTv ratio as a function of the annotation
#

gt = c[c$numVariants>minBinCutoff & c$category==0,]
novel = c[c$numVariants>minBinCutoff & c$category==1,]
dbsnp = c[c$numVariants>minBinCutoff & c$category==2,]

d = c[c$numVariants>minBinCutoff,]
cmin = min(d$titv)
cmax = max(d$titv)

outfile = paste(outputDir, "binnedTiTv.", annotationName, ".pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
plot(gt$value,gt$titv,xlab=annotationName,ylab="Ti/Tv ratio",pch=20);
m = weighted.mean(gt$value,gt$numVariants/sum(gt$numVariants))
ma = gt[gt$value > m,]
mb = gt[gt$value < m,]
m75 = weighted.mean(ma$value,ma$numVariants/sum(ma$numVariants))
m25 = weighted.mean(mb$value,mb$numVariants/sum(mb$numVariants))
abline(v=m,lty=2)
abline(v=m75,lty=2)
abline(v=m25,lty=2)
dev.off()

outfile = paste(outputDir, "binnedTiTv_novel.", annotationName, ".pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
plot(gt$value,gt$titv,xlab=annotationName,ylab="Ti/Tv ratio",pch=20,yim=c(cmin,cmax));
m = weighted.mean(gt$value,gt$numVariants/sum(gt$numVariants))
ma = gt[gt$value > m,]
mb = gt[gt$value < m,]
m75 = weighted.mean(ma$value,ma$numVariants/sum(ma$numVariants))
m25 = weighted.mean(mb$value,mb$numVariants/sum(mb$numVariants))
abline(v=m,lty=2)
abline(v=m75,lty=2)
abline(v=m25,lty=2)
points(novel$value,novel$titv,col="green",pch=20)
points(dbsnp$value,dbsnp$titv,col="blue",pch=20)
legend("topleft", c("overall","novel","dbsnp"),col=c("black","green","blue"),pch=c(20,20,20))
dev.off()


outfile = paste(outputDir, "binnedTiTv_quartiles.", annotationName, ".pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
plot(gt$value,gt$titv,xlab=annotationName,ylab="Ti/Tv ratio",pch=20,xlim=c(0,80))
abline(v=m,lty=2)
abline(v=m75,lty=2)
abline(v=m25,lty=2)
dev.off()

outfile = paste(outputDir, "binnedTiTv_hist.", annotationName, ".pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
plot(gt$value,gt$numVariants,xlab=annotationName,ylab="num Variants in bin",type="h");
dev.off()