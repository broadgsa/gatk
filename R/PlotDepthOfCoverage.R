args <- commandArgs(TRUE)
docBase <- args[1]

## APPEND THE SUFFIXES ##

locusStats <- paste(docBase,".sample_locus_statistics",sep="")
targetStats <- paste(docBase,".sample_interval_statistics",sep="")
sampleSum <- paste(docBase,".sample_summary_statistics",sep="")
sampleStats <- paste(docBase,".sample_statistics",sep="")
targetSum <- paste(docBase,".sample_interval_summary",sep="")

## DEFINE THE PLOTTING FUNCTIONS ##

PlotDepths <- function(X) {
	pdf("Depth_Histogram_All_Samples.pdf")
	Y <- as.matrix(X)
	colors <- rainbow(nrow(Y),gamma=0.8)
	plot(Y[1,],col=colors[1],type="b",xlab="",xaxt="n",ylab="Number of Loci")
	axis(1,labels=FALSE)
	labels <- colnames(X)
	text(1:ncol(Y),par("usr")[3]-(100/6000)*par("usr")[4],srt=45,adj=1,labels=labels,xpd=TRUE,cex=0.7)
	for ( jj in 2:nrow(Y) ) {
		points(Y[jj,],col=colors[jj],type="b")
	}
	ymax = par("usr")[4]
	xmax = par("usr")[2]
	legend(y=0.95*ymax,x=0.8*xmax,col=colors,rownames(X),lty=c(1),cex=0.5)
	dev.off()
}

PlotLocusQuantiles <- function(X) {
	pdf("Per_Sample_Coverage_Quantiles.pdf")
	Y <- as.matrix(X)
	Y <- Y/sum(Y[1,])
	Z <- matrix(nrow=nrow(Y),ncol=ncol(Y))
	for ( ii in 1:nrow(Y) ) {
		for ( jj in 1:ncol(Y) ) {
		# see how much density is in the remaining columns
			Z[ii,jj] = sum(Y[ii,jj:ncol(Y)])
		}	
	}
	
	medians = matrix(nrow=1,ncol=ncol(Z))
	quan90 = matrix(nrow=1,ncol=ncol(Z))
	for ( cc in 1:ncol(Z) ) {
		medians[cc] = median(Z[,cc])
		quan90[cc] = quantile(Z[,cc],0.9)
	}
	
	plot(t(medians),xlab="",xaxt="n",ylab="Proportion of loci with >X coverage",type="b",col="blue")
	axis(1,labels=FALSE)
	parseColNames <- function(K) {
		M = matrix(nrow=1,ncol=length(K))
		number = 0
		for ( lab in K ) {
			number = 1 + number
			g = unlist(strsplit(lab,split="_"))
			M[1,number] = g[2]
		}
		
		return(M)
	}
	labels <- parseColNames(colnames(X))
	text(1:length(labels),par("usr")[3]-0.025,srt=90,adj=1,labels=labels,xpd=TRUE,cex=(0.8/32)*length(labels),lheight=(0.8/32)*length(labels))
	points(t(quan90),type="b",col="red")
	legend(x=floor(0.6*length(labels)),y=1,c("50% of samples","90% of samples"),col=c("red","blue"),lty=c(1,1))
	dev.off()
}

HistogramMedians <- function(X) {
	pdf("Per_Sample_Median_Histogram.pdf")
	hist(as.numeric(as.matrix(unlist(X[1:nrow(X)-1,5]))),floor(nrow(X)/2),xlab="Median Coverage",ylab="Number of Samples", main="Median coverage acrosss samples",col="grey")
	dev.off()
}

HeatmapLocusTable <- function(X) {
	pdf("Locus_Coverage_HeatMap.pdf")
	Y <- as.matrix(X)
	heatmap(Y,Rowv=NA,Colv=NA)
	dev.off()
}

PlotMeanMedianQuartiles <- function(X) {
	pdf("Per_Sample_Mean_Quantile_Coverage.pdf")
	colors <- rainbow(4,start=0.6,end=0.9,gamma=1)
	means = as.numeric(as.matrix(unlist(X[1:nrow(X)-1,3])))
	medians = as.numeric(as.matrix(unlist(X[1:nrow(X)-1,5])))
	thirdQ = as.numeric(as.matrix(unlist(X[1:nrow(X)-1,4])))
	firstQ = as.numeric(as.matrix(unlist(X[1:nrow(X)-1,6])))
	plot(means,xlab="",ylab="Depth of Coverage",xaxt="n",col=colors[1],pch=3,type="b",ylim=c(0,max(thirdQ)))
	points(firstQ,col=colors[2],pch=2,type="b")
	points(medians,col=colors[3],pch=1,type="b")
	points(thirdQ,col=colors[4],pch=2,type="b")
	axis(1,labels=FALSE)
	labels <- X[1:nrow(X)-1,1]
	text(1:nrow(X)-1,par("usr")[3]-(50/2500)*par("usr")[4],srt=90,adj=1,labels=labels,xpd=TRUE,cex=0.5)
	text(5*nrow(X)/8,par("usr")[3]-(350/2500)*par("usr")[4],adj=1,labels="SAMPLE_ID",xpd=TRUE)
	legend(x=nrow(X)/10,y=par("usr")[4]-(200/2500)*par("usr")[4],c("Mean","25% Quantile","Median","75% Quantile"),col=colors,lty=c(1),cex=0.8,pch=c(3,2,1,2))
	dev.off()
}

PlotOnlyMeanMedian <- function(X) {
	pdf("Per_Sample_Mean_Median_Only.pdf")
	colors <- rainbow(2,start=0.6,end=0.9,gamma=1)
	means = as.numeric(as.matrix(unlist(X[1:nrow(X)-1,3])))
	medians = as.numeric(as.matrix(unlist(X[1:nrow(X)-1,5])))
	plot(means,xlab="",ylab="Depth of Coverage",xaxt="n",col=colors[1],pch=3,type="b",ylim=c(0,max(c(max(means),max(medians)))))
	points(medians,col=colors[2],pch=1,type="b")
	axis(1,labels=FALSE)
	labels <- X[1:nrow(X)-1,1]
	text(1:nrow(X)-1,par("usr")[3]-(50/2500)*par("usr")[4],srt=90,adj=1,labels=labels,xpd=TRUE,cex=0.5)
	text(5*nrow(X)/8,par("usr")[3]-(350/2500)*par("usr")[4],adj=1,labels="SAMPLE_ID",xpd=TRUE)
	legend(x=nrow(X)/10,y=par("usr")[4]-(200/2500)*par("usr")[4],c("Mean","Median"),col=colors,lty=c(1),cex=0.8,pch=c(3,2))
	dev.off()
}

PlotOnlyQuartiles <- function(X) {
	pdf("Per_Sample_Quartiles_Only.pdf")
	colors <- rainbow(2,start=0.6,end=0.9,gamma=1)
	thirdQ = as.numeric(as.matrix(unlist(X[1:nrow(X)-1,4])))
	firstQ = as.numeric(as.matrix(unlist(X[1:nrow(X)-1,6])))
	plot(thirdQ,xlab="",ylab="Depth of Coverage",xaxt="n",col=colors[1],pch=3,type="b",ylim=c(0,max(thirdQ)))
	points(firstQ,col=colors[2],pch=2,type="b")
	axis(1,labels=FALSE)
	labels <- X[1:nrow(X)-1,1]
	text(1:nrow(X)-1,par("usr")[3]-(50/2500)*par("usr")[4],srt=90,adj=1,labels=labels,xpd=TRUE,cex=0.5)
	text(5*nrow(X)/8,par("usr")[3]-(350/2500)*par("usr")[4],adj=1,labels="SAMPLE_ID",xpd=TRUE)
	legend(x=nrow(X)/10,y=par("usr")[4]-(200/2500)*par("usr")[4],c("75% Quantile","25% Quantile"),col=colors,lty=c(1),cex=0.8,pch=c(3,2))
	dev.off()
}

## PLOT SAMPLE STATISTICS
TO_PLOT <- read.table(sampleStats)
PlotDepths(TO_PLOT)
PlotLocusQuantiles(TO_PLOT)
## PLOT SAMPLE SUMMARY
TO_PLOT <- read.table(sampleSum,header=TRUE)
PlotMeanMedianQuartiles(TO_PLOT)
PlotOnlyMeanMedian(TO_PLOT)
PlotOnlyQuartiles(TO_PLOT)
HistogramMedians(TO_PLOT)
## PLOT LOCUS STATISTICS
TO_PLOT <- read.table(locusStats)
HeatmapLocusTable(TO_PLOT)
