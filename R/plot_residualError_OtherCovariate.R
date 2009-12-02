#!/broad/tools/apps/R-2.6.0/bin/Rscript

args <- commandArgs(TRUE)
verbose = TRUE

input = args[1]
covariateName = args[2]

#X11(width=7, height=14)
#outfile = paste(input, ".qual_diff_v_cycle.png", sep="")
#png(outfile, height=7, width=7, units="in", res=72) #height=1000, width=680)
outfile = paste(input, ".qual_diff_v_", covariateName, ".pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
c <- read.table(input, header=T)
c <- c[sort.list(c[,1]),]

d.good <- c[c$nBases >= 1000,]
d.1000 <- c[c$nBases < 1000,]
rmseBlue = sqrt(sum((d.good$Qempirical-d.good$Qreported)^2 * d.good$nBases) / sum(d.good$nBases) )
rmseAll = sqrt(sum((c$Qempirical-c$Qreported)^2 * c$nBases) / sum(c$nBases) )
theTitle = paste("RMSE_good = ", round(rmseBlue,digits=3), ", RMSE_all = ", round(rmseAll,digits=3))
if( length(d.good$nBases) == length(c$nBases) ) {
	theTitle = paste("RMSE = ", round(rmseAll,digits=3))
	}
if( is.numeric(c$Covariate) ) {
	plot(d.good$Covariate, d.good$Qempirical-d.good$Qreported, type="p", main=theTitle, ylab="Empirical - Reported Quality", 	xlab=covariateName, col="blue", pch=16, ylim=c(-10, 10), xlim=c(min(c$Covariate),max(c$Covariate)))
	points(d.1000$Covariate, d.1000$Qempirical-d.1000$Qreported, type="p", col="cornflowerblue", pch=16)
} else {
	plot(c$Covariate, c$Qempirical-c$Qreported, type="l", main=theTitle, ylab="Empirical - Reported Quality", 	xlab=covariateName, col="blue", ylim=c(-10, 10))
	points(d.1000$Covariate, d.1000$Qempirical-d.1000$Qreported, type="l", col="cornflowerblue")
}
dev.off()
#points(d.1000$Cycle, d.1000$Qempirical_Qreported, type="p", col="cornflowerblue", pch=16)

#
# Plot histogram of the covariate
#
e = d.good
f = d.1000
outfile = paste(input, ".", covariateName,"_hist.pdf", sep="")
pdf(outfile, height=7, width=7)
hst=subset(data.frame(e$Covariate, e$nBases), e.nBases != 0)
hst2=subset(data.frame(f$Covariate, f$nBases), f.nBases != 0)
if( is.numeric(c$Covariate) ) {
	plot(hst$e.Covariate, hst$e.nBases, type="h", lwd=4, main=paste(covariateName,"histogram"), xlab=covariateName, ylab="Count",yaxt="n")
	points(hst2$f.Covariate, hst2$f.nBases, type="h", lwd=4, col="cornflowerblue")
	axis(2,axTicks(2), format(axTicks(2), scientific=F))
} else {
	hst=subset(data.frame(c$Covariate, c$nBases), c.nBases != 0)
	plot(1:length(hst$c.Covariate), hst$c.nBases, type="h", lwd=4, main=paste(covariateName,"histogram"), xlab=covariateName, ylab="Count",yaxt="n",xaxt="n")
	axis(1, at=seq(1,length(hst$c.Covariate),2), labels = hst$c.Covariate[seq(1,length(hst$c.Covariate),2)])
	axis(2,axTicks(2), format(axTicks(2), scientific=F))
}
dev.off()