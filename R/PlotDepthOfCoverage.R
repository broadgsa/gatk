#!/broad/tools/apps/R-2.6.0/bin/Rscript

args <- commandArgs(TRUE)
docBase <- args[1]
locusStats <- cat(docBase,".sample_locus_statistics")
targetStats <- cat(docBase,".sample_interval_statistics")
sampleSum <- cat(docBase,".sample_summary_statistics")
sampleStats <- cat(docBase,".sample_statistics")
targetSum <- cat(docBase(".sample_interval_summary"))

PlotDepths <- function(X) {
	Y <- as.matrix(X)
	colors <- rainbow(nrow(Y),gamma=0.8)
	plot(Y[1,],col=colors[1],type="b",xlab="",xaxt="n",ylab="Number of Loci")
	axis(1,labels=FALSE)
	labels <- colnames(X)
	text(1:ncol(Y),par("usr")[3]-0.25,srt=45,adj=1,labels=labels,xpd=TRUE)
	for ( jj in 2:nrow(Y) ) {
		points(Y[jj,],col=colors[jj],type="b")
	}
	ymax = par("usr")[4]
	xmax = par("usr")[2]
	legend(y=0.95*ymax,x=0.8*xmax,col=colors,rownames(X),lty=c(1))
}

PlotLocusQuantiles <- function(X) {
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
		M = matrix(nrow=1,ncol=length(labels))
		number = 0
		for ( lab in K ) {
			number = 1 + number
			M[1,number] = unlist(strsplit(lab,split="_"))[2]
		}
		
		return(M)
	}
	labels <- parseColNames(colnames(X))
	text(1:length(labels),par("usr")[3]-0.05,srt=90,adj=1,labels=labels,xpd=TRUE)
	points(t(quan90),type="b",col="red")
	legend(x=floor(0.6*length(labels)),y=1,c("50% of samples","90% of samples"),col=c("red","blue"),lty=c(1,1))
}

HistogramMedians <- function(X) {
	hist(as.numeric(as.matrix(unlist(X[1:nrow(X)-1,5]))),floor(nrow(X)/2),xlab="Median Coverage",ylab="Number of Samples", main="Median coverage acrosss samples",col="grey")
}

HeatmapLocusTable <- function(X) {
	Y <- as.matrix(X)
	heatmap(Y,Rowv=NA,Colv=NA)
}

PlotMeanMeadianQuartiles <- function(X) {
	colors <- rainbow(4,start=0.6,end=0.9,gamma=1)
	means = as.numeric(as.matrix(unlist(X[1:nrow(X)-1,3])))
	medians = as.numeric(as.matrix(unlist(X[1:nrow(X)-1,5])))
	thirdQ = as.numeric(as.matrix(unlist(X[1:nrow(X)-1,4])))
	firstQ = as.numeric(as.matrix(unlist(X[1:nrow(X)-1,6])))
	plot(means,xlab="",ylab="Depth of Coverage",xaxt="n",col=colors[1],pch=3,type="b")
	points(firstQ,col=colors[2],pch=2,type="b")
	points(medians,col=colors[3],pch=1,type="b")
	points(thirdQ,col=colors[4],pch=2,type="b")
	axis(1,labels=FALSE)
	labels <- X[1:nrow(X)-1,1]
	text(1:nrow(X)-1,par("usr")[3]-0.25,srt=45,adj=1,labels=labels,xpd=TRUE)
}