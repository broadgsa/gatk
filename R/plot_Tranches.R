#!/bin/env Rscript

args <- commandArgs(TRUE)
verbose = TRUE

input = args[1]
targetTITV = as.numeric(args[2])

# -----------------------------------------------------------------------------------------------
# Useful general routines
# -----------------------------------------------------------------------------------------------

MIN_FP_RATE = 0.01

titvFPEst <- function(titvExpected, titvObserved) { 
    max(min(1 - (titvObserved - 0.5) / (titvExpected - 0.5), 1), MIN_FP_RATE) 
}

titvFPEstV <- function(titvExpected, titvs) {
    sapply(titvs, function(x) titvFPEst(titvExpected, x))
}

nTPFP <- function(nVariants, FDR) {
    return(list(TP = nVariants * (1 - FDR/100), FP = nVariants * (FDR / 100)))
}

leftShift <- function(x, leftValue = 0) {
    r = rep(leftValue, length(x))
    for ( i in 1:(length(x)-1) ) {
        #print(list(i=i))
        r[i] = x[i+1]
    }
    r
}

# -----------------------------------------------------------------------------------------------
# Tranches plot
# -----------------------------------------------------------------------------------------------
data2 = read.table(paste(input,".tranches",sep=""),sep=",",head=T)
cols = c("cornflowerblue", "cornflowerblue", "darkorange", "darkorange")
density=c(20, -1, -1, 20)
outfile = paste(input, ".FDRtranches.pdf", sep="")
pdf(outfile, height=7, width=8)
alpha = 1 - titvFPEstV(targetTITV, data2$novelTITV)
#print(alpha)

numGood = round(alpha * data2$numNovel);

#numGood = round(data2$numNovel * (1-data2$FDRtranche/100))
numBad = data2$numNovel - numGood;

numPrevGood = leftShift(numGood, 0)
numNewGood = numGood - numPrevGood
numPrevBad = leftShift(numBad, 0)
numNewBad = numBad - numPrevBad

d=matrix(c(numPrevGood,numNewGood, numNewBad, numPrevBad),4,byrow=TRUE)
#print(d)
barplot(d/1000,horiz=TRUE,col=cols,space=0.2,xlab="Number of Novel Variants (1000s)",ylab="Novel Ti/Tv   -->   FDR (%)", density=density) # , xlim=c(250000,350000))
#abline(v= d[2,dim(d)[2]], lty=2)
#abline(v= d[1,3], lty=2)
legend(10000/1000, 2.25, c('Cumulative TPs','Tranch-specific TPs', 'Tranch-specific FPs', 'Cumulative FPs' ), fill=cols, density=density, bg='white', cex=1.25)
axis(2,line=-1,at=0.7+(0:(length(data2$FDRtranche)-1))*1.2,tick=FALSE,labels=data2$FDRtranche)
axis(2,line=0.4,at=0.7+(0:(length(data2$FDRtranche)-1))*1.2,tick=FALSE,labels=data2$novelTITV)
dev.off()


#
#data2 = read.table(paste(input,".tranches",sep=""),sep=",",head=T)
#cols = c("steelblue","orange")
#outfile = paste(input, ".FDRtranches.pdf", sep="")
#pdf(outfile, height=7, width=8)
#alpha = (data2$novelTITV - 0.5) / (targetTITV - 0.5);
#numGood = round(alpha * data2$numNovel);
#numBad = data2$numNovel - numGood;
#d=matrix(c(numGood,numBad),2,byrow=TRUE)
#barplot(d,horiz=TRUE,col=cols,space=0.2,xlab="Number of Novel Variants",ylab="Novel Ti/Tv   -->   FDR (%)")
#legend('topright',c('implied TP','implied FP'),col=cols,lty=1,lwd=16)
#axis(2,line=-1,at=0.7+(0:(length(data2$FDRtranche)-1))*1.2,tick=FALSE,labels=data2$FDRtranche)
#axis(2,line=0.4,at=0.7+(0:(length(data2$FDRtranche)-1))*1.2,tick=FALSE,labels=data2$novelTITV)
#dev.off()
