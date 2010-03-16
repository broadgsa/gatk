#!/broad/tools/apps/R-2.6.0/bin/Rscript

args <- commandArgs(TRUE)
verbose = TRUE

input = args[1]
targetTITV = as.numeric(args[2])


data = read.table(input,sep=",",head=T)
maxVars = max(data$numKnown, data$numNovel)
maxTITV = max(data$knownTITV[is.finite(data$knownTITV) & data$numKnown>2000], data$novelTITV[is.finite(data$novelTITV) & data$numNovel > 2000], targetTITV)
minTITV = min(data$knownTITV[length(data$knownTITV)], data$novelTITV[length(data$novelTITV)], targetTITV)


outfile = paste(input, ".optimizationCurve.pdf", sep="")
pdf(outfile, height=7, width=8)

par(mar=c(4,4,1,4),cex=1.3)
plot(data$pCut, data$knownTITV, axes=F,xlab="Keep variants with p(true) > X",ylab="",ylim=c(minTITV,maxTITV),col="Blue",pch=20)
points(data$pCut, data$novelTITV,,col="DarkBlue",pch=20)
abline(h=targetTITV,lty=3,col="Blue")
axis(side=2,col="DarkBlue")
axis(side=1)
mtext("Ti/Tv Ratio", side=2, line=2, col="blue",cex=1.4)
legend("left", c("Known Ti/Tv","Novel Ti/Tv"), col=c("Blue","DarkBlue"), pch=c(20,20),cex=0.7)
par(new=T)
plot(data$pCut, data$numKnown, axes=F,xlab="",ylab="",ylim=c(0,maxVars),col="Green",pch=20)
points(data$pCut, data$numNovel,col="DarkGreen",pch=20)
axis(side=4,col="DarkGreen")
mtext("Number of Variants", side=4, line=2, col="DarkGreen",cex=1.4)
legend("topright", c("Known","Novel"), col=c("Green","DarkGreen"), pch=c(20,20),cex=0.7)
dev.off()