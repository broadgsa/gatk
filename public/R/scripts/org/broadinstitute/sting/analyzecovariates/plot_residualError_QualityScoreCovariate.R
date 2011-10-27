#!/bin/env Rscript

library(tools)

args <- commandArgs(TRUE)

input = args[1]
Qcutoff = as.numeric(args[2])
maxQ = as.numeric(args[3])
maxHist = as.numeric(args[4])

t=read.table(input, header=T)

#
# Plot of reported quality versus empirical quality
#

outfile = paste(input, ".quality_emp_v_stated.pdf", sep="")
pdf(outfile, height=7, width=7)
d.good <- t[t$nBases >= 10000 & t$Qreported >= Qcutoff,]
d.1000 <- t[t$nBases < 1000  & t$Qreported >= Qcutoff,]
d.10000 <- t[t$nBases < 10000 & t$nBases >= 1000  & t$Qreported >= Qcutoff,]
f <- t[t$Qreported < Qcutoff,]
e <- rbind(d.good, d.1000, d.10000)
rmseGood = sqrt( sum(as.numeric((d.good$Qempirical-d.good$Qreported)^2 * d.good$nBases)) / sum(as.numeric(d.good$nBases)) ) # prevent integer overflow with as.numeric, ugh
rmseAll = sqrt( sum(as.numeric((e$Qempirical-e$Qreported)^2 * e$nBases)) / sum(as.numeric(e$nBases)) )
theTitle = paste("RMSE_good =", round(rmseGood,digits=3), ", RMSE_all =", round(rmseAll,digits=3))
if( length(t$nBases) - length(f$nBases) == length(d.good$nBases) ) {
	theTitle = paste("RMSE =", round(rmseAll,digits=3));
}
plot(d.good$Qreported, d.good$Qempirical, type="p", col="blue", main=theTitle, xlim=c(0,maxQ), ylim=c(0,maxQ), pch=16, xlab="Reported quality score", ylab="Empirical quality score")
points(d.1000$Qreported, d.1000$Qempirical, type="p", col="lightblue", pch=16)
points(d.10000$Qreported, d.10000$Qempirical, type="p", col="cornflowerblue", pch=16)
points(f$Qreported, f$Qempirical, type="p", col="maroon1", pch=16)
abline(0,1, lty=2)
dev.off()

if (exists('compactPDF')) {
  compactPDF(outfile)
}

#
# Plot Q empirical histogram
#

outfile = paste(input, ".quality_emp_hist.pdf", sep="")
pdf(outfile, height=7, width=7)
hst=subset(data.frame(e$Qempirical, e$nBases), e.nBases != 0)
hst2=subset(data.frame(f$Qempirical, f$nBases), f.nBases != 0)
percentBases=hst$e.nBases / sum(as.numeric(hst$e.nBases))
entropy = -sum(log2(percentBases)*percentBases)
yMax = max(hst$e.nBases)
if(maxHist != 0) {
yMax = maxHist
}
plot(hst$e.Qempirical, hst$e.nBases, type="h", lwd=4, xlim=c(0,maxQ), ylim=c(0,yMax), main=paste("Empirical quality score histogram, entropy = ",round(entropy,digits=3)), xlab="Empirical quality score", ylab="Number of Bases",yaxt="n")
points(hst2$f.Qempirical, hst2$f.nBases, type="h", lwd=4, col="maroon1")
axis(2,axTicks(2), format(axTicks(2), scientific=F))
dev.off()

if (exists('compactPDF')) {
  compactPDF(outfile)
}

#
# Plot Q reported histogram
#

outfile = paste(input, ".quality_rep_hist.pdf", sep="")
pdf(outfile, height=7, width=7)
hst=subset(data.frame(e$Qreported, e$nBases), e.nBases != 0)
hst2=subset(data.frame(f$Qreported, f$nBases), f.nBases != 0)
yMax = max(hst$e.nBases)
if(maxHist != 0) {
yMax = maxHist
}
plot(hst$e.Qreported, hst$e.nBases, type="h", lwd=4, xlim=c(0,maxQ), ylim=c(0,yMax), main=paste("Reported quality score histogram, entropy = ",round(entropy,digits=3)), xlab="Reported quality score", ylab="Number of Bases",yaxt="n")
points(hst2$f.Qreported, hst2$f.nBases, type="h", lwd=4, col="maroon1")
axis(2,axTicks(2), format(axTicks(2), scientific=F))
dev.off()

if (exists('compactPDF')) {
  compactPDF(outfile)
}
