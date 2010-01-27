#!/broad/tools/apps/R-2.6.0/bin/Rscript

args <- commandArgs(TRUE)
verbose = TRUE

input = args[1]
annotationName = args[2]
minBinCutoff = as.numeric(args[3])


c <- read.table(input, header=T)

#
# Plot TiTv ratio as a function of the annotation
#

all = c[c$numVariants>minBinCutoff & c$category=="all",]
novel = c[c$numVariants>minBinCutoff & c$category=="novel",]
dbsnp = c[c$numVariants>minBinCutoff & c$category=="dbsnp",]

d = c[c$numVariants>minBinCutoff,]
ymin = min(d$titv)
ymax = max(d$titv)
xmin = min(d$value)
xmax = max(d$value)

outfile = paste(input, ".TiTv.", annotationName, ".pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
plot(all$value,all$titv,xlab=annotationName,ylab="Ti/Tv ratio",pch=20,ylim=c(ymin,ymax),xaxt="n",ps=14);
axis(1,axTicks(1), format(axTicks(1), scientific=F))
m = weighted.mean(all$value,all$numVariants/sum(all$numVariants))
ma = all[all$value > m,]
mb = all[all$value < m,]
m75 = weighted.mean(ma$value,ma$numVariants/sum(ma$numVariants))
m25 = weighted.mean(mb$value,mb$numVariants/sum(mb$numVariants))
abline(v=m,lty=2)
abline(v=m75,lty=2)
abline(v=m25,lty=2)
points(novel$value,novel$titv,col="green",pch=20)
points(dbsnp$value,dbsnp$titv,col="blue",pch=20)
legend("topleft", c("all","novel","dbsnp"),col=c("black","green","blue"),pch=c(20,20,20))
dev.off()

outfile = paste(input, ".TiTv_log.", annotationName, ".pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
plot(all$value,all$titv,xlab=annotationName,log="x",ylab="Ti/Tv ratio",pch=20,ylim=c(ymin,ymax),xaxt="n",ps=14);
axis(1,axTicks(1), format(axTicks(1), scientific=F))
abline(v=m,lty=2)
abline(v=m75,lty=2)
abline(v=m25,lty=2)
points(novel$value,novel$titv,col="green",pch=20)
points(dbsnp$value,dbsnp$titv,col="blue",pch=20)
legend("topleft", c("all","novel","dbsnp"),col=c("black","green","blue"),pch=c(20,20,20))
dev.off()


outfile = paste(input, "TiTv_hist.", annotationName, ".pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
plot(all$value,all$numVariants,xlab=annotationName,ylab="num variants in bin",type="h",xaxt="n",ps=14);
axis(1,axTicks(1), format(axTicks(1), scientific=F))
dev.off()