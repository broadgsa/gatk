#!/broad/tools/apps/R-2.6.0/bin/Rscript

args <- commandArgs(TRUE)
verbose = TRUE

input = args[1]
annotationName = args[2]
minBinCutoff = as.numeric(args[3])
medianNumVariants = args[4]

c <- read.table(input, header=T)

all = c[c$numVariants>minBinCutoff & c$category=="all",]
novel = c[c$numVariants>minBinCutoff & c$category=="novel",]
dbsnp = c[c$numVariants>minBinCutoff & c$category=="dbsnp",]
truth = c[c$numVariants>minBinCutoff & c$category=="truth",]

#
# Calculate min, max, medians
#

d = c[c$numVariants>minBinCutoff,]
ymin = min(d$titv)
ymax = max(d$titv)
xmin = min(d$value)
xmax = max(d$value)
m = weighted.mean(all$value,all$numVariants/sum(all$numVariants))
ma = all[all$value > m,]
mb = all[all$value < m,]
m75 = weighted.mean(ma$value,ma$numVariants/sum(ma$numVariants))
m25 = weighted.mean(mb$value,mb$numVariants/sum(mb$numVariants))
if(medianNumVariants == "true") {
vc = cumsum( all$numVariants/sum(all$numVariants) )
m25 = all$value[ max(which(vc<=0.25)) ]
m = all$value[ max(which(vc<=0.5)) ]
m75 = all$value[ min(which(vc>=0.75)) ]
}

#
# Plot TiTv ratio as a function of the annotation
#

outfile = paste(input, ".TiTv.pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
plot(all$value,all$titv,xlab=annotationName,ylab="Ti/Tv Ratio",pch=20,ylim=c(ymin,ymax),xaxt="n",ps=14);
axis(1,axTicks(1), format(axTicks(1), scientific=F))
abline(v=m,lty=2)
abline(v=m75,lty=2)
abline(v=m25,lty=2)
points(novel$value,novel$titv,col="green",pch=20)
points(dbsnp$value,dbsnp$titv,col="blue",pch=20)
if( sum(all$truePositive==0) != length(all$truePositive) ) {
points(truth$value,truth$titv,col="magenta",pch=20)
legend("topleft", c("all","novel","dbsnp","truth"),col=c("black","green","blue","magenta"),pch=c(20,20,20,20))
} else {
legend("topleft", c("all","novel","dbsnp"),col=c("black","green","blue"),pch=c(20,20,20))
}
dev.off()

#
# Plot TiTv ratio as a function of the annotation, log scale on the x-axis
#

outfile = paste(input, ".TiTv_log.pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
plot(all$value,all$titv,xlab=annotationName,log="x",ylab="Ti/Tv Ratio",pch=20,ylim=c(ymin,ymax),xaxt="n",ps=14);
axis(1,axTicks(1), format(axTicks(1), scientific=F))
abline(v=m,lty=2)
abline(v=m75,lty=2)
abline(v=m25,lty=2)
points(novel$value,novel$titv,col="green",pch=20)
points(dbsnp$value,dbsnp$titv,col="blue",pch=20)
if( sum(all$truePositive==0) != length(all$truePositive) ) {
points(truth$value,truth$titv,col="magenta",pch=20)
legend("topleft", c("all","novel","dbsnp","truth"),col=c("black","green","blue","magenta"),pch=c(20,20,20,20))
} else {
legend("topleft", c("all","novel","dbsnp"),col=c("black","green","blue"),pch=c(20,20,20))
}
dev.off()

#
# Plot dbsnp and true positive rate as a function of the annotation
#

ymin = min(all$dbsnp)
ymax = max(all$dbsnp)
outfile = paste(input, ".truthRate.pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
yLabel = "DBsnp Rate"
if( sum(all$truePositive==0) != length(all$truePositive) ) {
t = all[all$truePositive>0,]
yLabel = "DBsnp/True Positive Rate"
ymin = min(min(all$dbsnp),min(t$truePositive))
ymax = max(max(all$dbsnp),max(t$truePositive))
}
plot(all$value,all$dbsnp,xlab=annotationName,ylab=yLabel,pch=20,ylim=c(ymin,ymax),xaxt="n",ps=14);
axis(1,axTicks(1), format(axTicks(1), scientific=F))
abline(v=m,lty=2)
abline(v=m75,lty=2)
abline(v=m25,lty=2)
if( sum(all$truePositive==0) != length(all$truePositive) ) {
points(t$value,t$truePositive,col="magenta",pch=20);
legend("topleft", c("dbsnp","truth"),col=c("black","magenta"),pch=c(20,20))
}
dev.off()

#
# Plot dbsnp and true positive rate as a function of the annotation, log scale on the x-axis
#

outfile = paste(input, ".truthRate_log.pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
yLabel = "DBsnp Rate"
if( sum(all$truePositive==0) != length(all$truePositive) ) {
yLabel = "DBsnp/Truth Rate"
}
plot(all$value,all$dbsnp,xlab=annotationName,log="x",ylab=yLabel,ylim=c(ymin,ymax),pch=20,xaxt="n",ps=14);
axis(1,axTicks(1), format(axTicks(1), scientific=F))
abline(v=m,lty=2)
abline(v=m75,lty=2)
abline(v=m25,lty=2)
if( sum(all$truePositive==0) != length(all$truePositive) ) {
points(t$value,t$truePositive,col="magenta",pch=20);
legend("topleft", c("dbsnp","truth"),col=c("black","magenta"),pch=c(20,20))
}
dev.off()

#
# Plot histogram of the annotation's value
#

outfile = paste(input, ".Histogram.pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
plot(all$value,all$numVariants,xlab=annotationName,ylab="Num variants in bin",type="h",xaxt="n",ps=14,lwd=4);
axis(1,axTicks(1), format(axTicks(1), scientific=F))
dev.off()

#
# Plot histogram of the annotation's value, log scale on x-axis
#

outfile = paste(input, ".Histogram_log.pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
plot(all$value,all$numVariants,xlab=annotationName,log="x",ylab="Num variants in bin",type="h",xaxt="n",ps=14,lwd=4);
axis(1,axTicks(1), format(axTicks(1), scientific=F))
dev.off()