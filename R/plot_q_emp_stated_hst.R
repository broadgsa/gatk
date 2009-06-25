#!/broad/tools/apps/R-2.6.0/bin/Rscript

args <- commandArgs(TRUE)

input = args[1]

t=read.table(input, header=T)
#t=read.csv(input)
#par(mfrow=c(2,1), cex=1.2)

#outfile = paste(input, ".quality_emp_v_stated.png", sep="")
#png(outfile, height=7, width=7, units="in", res=72) # height=1000, width=446)
outfile = paste(input, ".quality_emp_v_stated.pdf", sep="")
pdf(outfile, height=7, width=7)
d.good <- t[t$nMismatches >= 1000,]
d.100 <- t[t$nMismatches < 100,]
d.1000 <- t[t$nMismatches < 1000 & t$nMismatches >= 100,]
plot(d.good$Qreported, d.good$Qempirical, type="p", col="blue", xlim=c(0,45), ylim=c(0,45), pch=16, xlab="Reported quality score", ylab="Empirical quality score", main="Reported vs. empirical quality scores")
points(d.100$Qreported, d.100$Qempirical, type="p", col="lightblue", pch=16)
points(d.1000$Qreported, d.1000$Qempirical, type="p", col="cornflowerblue", pch=16)
abline(0,1, lty=2)
dev.off()

#outfile = paste(input, ".quality_emp_hist.png", sep="")
#png(outfile, height=7, width=7, units="in", res=72) # height=1000, width=446)
outfile = paste(input, ".quality_emp_hist.pdf", sep="")
pdf(outfile, height=7, width=7)
hst=subset(data.frame(t$Qempirical, t$nBases), t.nBases != 0)
plot(hst$t.Qempirical, hst$t.nBases, type="h", lwd=3, xlim=c(0,45), main="Reported quality score histogram", xlab="Empirical quality score", ylab="Count", yaxt="n")
axis(2,axTicks(2), format(axTicks(2), scientific=F))
dev.off()
