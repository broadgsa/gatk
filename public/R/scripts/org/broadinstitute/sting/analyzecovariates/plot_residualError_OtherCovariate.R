#!/bin/env Rscript

library(tools)

args <- commandArgs(TRUE)
verbose = TRUE

input = args[1]
covariateName = args[2]

outfile = paste(input, ".qual_diff_v_", covariateName, ".pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
c <- read.table(input, header=T)
c <- c[sort.list(c[,1]),]

#
# Plot residual error as a function of the covariate
#

d.good <- c[c$nBases >= 1000,]
d.1000 <- c[c$nBases < 1000,]
rmseGood = sqrt( sum(as.numeric((d.good$Qempirical-d.good$Qreported)^2 * d.good$nBases)) / sum(as.numeric(d.good$nBases)) ) # prevent integer overflow with as.numeric, ugh
rmseAll = sqrt( sum(as.numeric((c$Qempirical-c$Qreported)^2 * c$nBases)) / sum(as.numeric(c$nBases)) )
theTitle = paste("RMSE_good =", round(rmseGood,digits=3), ", RMSE_all =", round(rmseAll,digits=3))
if( length(d.good$nBases) == length(c$nBases) ) {
	theTitle = paste("RMSE =", round(rmseAll,digits=3))
}
# Don't let residual error go off the edge of the plot
d.good$residualError = d.good$Qempirical-d.good$Qreported
d.good$residualError[which(d.good$residualError > 10)] = 10
d.good$residualError[which(d.good$residualError < -10)] = -10
d.1000$residualError = d.1000$Qempirical-d.1000$Qreported
d.1000$residualError[which(d.1000$residualError > 10)] = 10
d.1000$residualError[which(d.1000$residualError < -10)] = -10
c$residualError = c$Qempirical-c$Qreported
c$residualError[which(c$residualError > 10)] = 10
c$residualError[which(c$residualError < -10)] = -10
pointType = "p"
if( length(c$Covariate) <= 20 ) {
    pointType = "o"
}
if( is.numeric(c$Covariate) ) {
	plot(d.good$Covariate, d.good$residualError, type=pointType, main=theTitle, ylab="Empirical - Reported Quality", xlab=covariateName, col="blue", pch=20, ylim=c(-10, 10), xlim=c(min(c$Covariate),max(c$Covariate)))
	points(d.1000$Covariate, d.1000$residualError, type=pointType, col="cornflowerblue", pch=20)
} else { # Dinuc (and other non-numeric covariates) are different to make their plots look nice
	plot(c$Covariate, c$residualError, type="l", main=theTitle, ylab="Empirical - Reported Quality", xlab=covariateName, col="blue", ylim=c(-10, 10))
	points(d.1000$Covariate, d.1000$residualError, type="l", col="cornflowerblue")
}
dev.off()

if (exists('compactPDF')) {
  compactPDF(outfile)
}

#
# Plot mean quality versus the covariate
#

outfile = paste(input, ".reported_qual_v_", covariateName, ".pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
pointType = "p"
if( length(c$Covariate) <= 20 ) {
    pointType = "o"
}
theTitle = paste("Quality By", covariateName);
if( is.numeric(c$Covariate) ) {
	plot(d.good$Covariate, d.good$Qreported, type=pointType, main=theTitle, ylab="Mean Reported Quality", xlab=covariateName, col="blue", pch=20, ylim=c(0, 40), xlim=c(min(c$Covariate),max(c$Covariate)))
	points(d.1000$Covariate, d.1000$Qreported, type=pointType, col="cornflowerblue", pch=20)
} else { # Dinuc (and other non-numeric covariates) are different to make their plots look nice
	plot(c$Covariate, c$Qreported, type="l", main=theTitle, ylab="Mean Reported Quality", xlab=covariateName, col="blue", ylim=c(0, 40))
	points(d.1000$Covariate, d.1000$Qreported, type="l", col="cornflowerblue")
}
dev.off()

if (exists('compactPDF')) {
  compactPDF(outfile)
}

#
# Plot histogram of the covariate
#

e = d.good
f = d.1000
outfile = paste(input, ".", covariateName,"_hist.pdf", sep="")
pdf(outfile, height=7, width=7)
hst=subset(data.frame(e$Covariate, e$nBases), e.nBases != 0)
hst2=subset(data.frame(f$Covariate, f$nBases), f.nBases != 0)

lwdSize=2
if( length(c$Covariate) <= 20 ) {
    lwdSize=7
} else if( length(c$Covariate) <= 70 ) {
    lwdSize=4
}

if( is.numeric(c$Covariate) ) {
    if( length(hst$e.Covariate) == 0 ) {
        plot(hst2$f.Covariate, hst2$f.nBases, type="h", lwd=lwdSize, col="cornflowerblue", main=paste(covariateName,"histogram"), ylim=c(0, max(hst2$f.nBases)), xlab=covariateName, ylab="Count",yaxt="n",xlim=c(min(c$Covariate),max(c$Covariate)))
    } else {
	    plot(hst$e.Covariate, hst$e.nBases, type="h", lwd=lwdSize, main=paste(covariateName,"histogram"), xlab=covariateName, ylim=c(0, max(hst$e.nBases)),ylab="Number of Bases",yaxt="n",xlim=c(min(c$Covariate),max(c$Covariate)))
	    points(hst2$f.Covariate, hst2$f.nBases, type="h", lwd=lwdSize, col="cornflowerblue")
	}
	axis(2,axTicks(2), format(axTicks(2), scientific=F))
} else { # Dinuc (and other non-numeric covariates) are different to make their plots look nice
	hst=subset(data.frame(c$Covariate, c$nBases), c.nBases != 0)
	plot(1:length(hst$c.Covariate), hst$c.nBases, type="h", lwd=lwdSize, main=paste(covariateName,"histogram"), ylim=c(0, max(hst$c.nBases)),xlab=covariateName, ylab="Number of Bases",yaxt="n",xaxt="n")
	if( length(hst$c.Covariate) > 9 ) {
	    axis(1, at=seq(1,length(hst$c.Covariate),2), labels = hst$c.Covariate[seq(1,length(hst$c.Covariate),2)])
	} else {
	    axis(1, at=seq(1,length(hst$c.Covariate),1), labels = hst$c.Covariate)
	}
	axis(2,axTicks(2), format(axTicks(2), scientific=F))
}
dev.off()

if (exists('compactPDF')) {
  compactPDF(outfile)
}
