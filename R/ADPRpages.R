#These functions each make a page for the ADPR. They assume a pdf with the following parameters for best formatting:
#pdf(file=paste(sample_sets, ".pdf", sep=""), width=22, height=15, pagecentre=TRUE, pointsize=24)        


library(gplots)
library(ReadImages)

##defaults<-par(no.readonly = TRUE)


tearsheet<-function(lanetable, sampletable, variant, Protocol, Sequencer){
	
	#define layout
	layout(matrix(c(1,1,2,4,3,5), ncol=2, nrow=3, byrow=TRUE), heights=c(1, 2.5,2.5,), respect=FALSE)
    
    #prep for title bar
    title=paste(sample_sets, ": TEAR SHEET", sep="")
    drop<-read.jpeg("tearsheetdrop.jpg")
    
    #plot title bar
    par(mar=c(0,0,0,0))
    plot(drop)
    text(100, 40, title, family="serif", adj=c(0,0), cex=3, col=gray(.25))
    

	#calc by lane stuff
	sdlane<-rep("NA", 6)
	meanlane<-sdlane
    
    attach(lanetable);
	   	
	callable.target<-HS_TARGET_TERRITORY[1]; 
    singlelanes<-length(which(Lane.Type=="Single"));
    pairedlanes<-length(which(Lane.Type=="Paired"));
   	meanlane[1]<-round(mean(AL_TOTAL_READS, na.rm=TRUE)/10^6, 2);
	sdlane[1]<-round(sd(AL_TOTAL_READS, na.rm=TRUE)/10^6, 2);
	meanlane[2]<-round(mean(HS_ON_TARGET_BASES, na.rm=TRUE)/10^6, 2);
	sdlane[2]<-round(sd(HS_ON_TARGET_BASES, na.rm=TRUE)/10^6, 2);
	meanlane[3]<-round(mean(HS_MEAN_TARGET_COVERAGE, na.rm=TRUE));
    sdlane[3]<-round(sd(HS_MEAN_TARGET_COVERAGE, na.rm=TRUE));
    meanlane[4]<-round(mean(HS_PCT_TARGET_BASES_10X, na.rm=TRUE));
    meanlane[5]<-round(mean(HS_PCT_TARGET_BASES_20X, na.rm=TRUE));
    meanlane[6]<-round(mean(HS_PCT_TARGET_BASES_30X, na.rm=TRUE));
    sdlane[4]<-round(sd(HS_PCT_TARGET_BASES_10X, na.rm=TRUE));
    sdlane[5]<-round(sd(HS_PCT_TARGET_BASES_20X, na.rm=TRUE));
    sdlane[6]<-round(sd(HS_PCT_TARGET_BASES_30X, na.rm=TRUE))
	   		
	names<-paste(Flowcell, "-", Lane, sep="")
       
    detach(lanetable)
			
	meansamp<-rep("NA", 6)
	sdsamp<-meansamp

	#Calc by sample metrics   	
    attach(bysample);
	baits<-Bait.Set[1]
	alllanes<-signif(sum(X..Lanes.included.in.aggregation, na.rm = TRUE))
	mean.lanes.samp<-signif(mean(X..Lanes.included.in.aggregation, na.rm = TRUE));
 	sd.lanes.samp<-signif(sd(X..Lanes.included.in.aggregation, na.rm=TRUE));
	mean.mrl.samp<-signif(mean(Mean.Read.Length, na.rm=TRUE));
	sd.mrl.samp<-signif(sd(Mean.Read.Length, na.rm=TRUE));
	meansamp[1]<-round(mean(Total.Reads, na.rm=TRUE)/10^6, 2);
	sdsamp[1]<-round(sd(Total.Reads, na.rm=TRUE)/10^6, 2);
	meansamp[2]<-round(mean(On.Target.Bases..HS., na.rm=TRUE)/10^6, 2);
	sdsamp[2]<-round(sd(On.Target.Bases..HS., na.rm=TRUE)/10^6, 2);
	meansamp[3]<-round(mean(Mean.Target.Coverage..HS., na.rm=TRUE));
	sdsamp[3]<-round(sd(Mean.Target.Coverage..HS., na.rm=TRUE));
	meansamp[4]<-round(mean(PCT.Target.Bases.10x..HS., na.rm=TRUE));
	meansamp[5]<-round(mean(PCT.Target.Bases.20x..HS., na.rm=TRUE));
	meansamp[6]<-round(mean(PCT.Target.Bases.30x..HS., na.rm=TRUE));
	sdsamp[4]<-round(sd(PCT.Target.Bases.10x..HS., na.rm=TRUE));
	sdsamp[5]<-round(sd(PCT.Target.Bases.20x..HS., na.rm=TRUE));
	sdsamp[6]<-round(sd(PCT.Target.Bases.30x..HS., na.rm=TRUE));
			
	detach(bysample);

	#calc variant stuff
	attach(variant)
	SNPS<-c(ti_count[which(filter_name=="called")]+tv_count[which(filter_name=="called")])
	titvs<-c(ti.tv_ratio[which(filter_name=="called")])
	detach(variant)
	
	#prep stuff. 
	summary<-c(nrow(bysample), Protocol, baits, paste(callable.target, "bases"))
	summary2<-c(Sequencer, alllanes, paste(mean.lanes.samp, "+/-", sd.lanes.samp), paste(singlelanes, "single lanes,", pairedlanes, "paired lanes"), paste(mean.mrl.samp, "+/-", sd.mrl.samp))
	samps<-paste(meansamp, c("M", "M", "x", "%", "%", "%"), " +/- ", sdsamp, c("M", "M", "x", "%", "%", "%"), sep="")
	lanes<-paste(meanlane, c("M", "M", "x", "%", "%", "%"), " +/- ", sdlane, c("M", "M", "x", "%", "%", "%"), sep="")
			
	#print out 4 tables in R
	table1<-cbind(summary)
	rownames(table1)<-c("Samples","Sequencing Protocol", "Bait Design","Callable Target")
	par(mar=c(4,4,4,4))
	textplot(table1, col.rownames="darkblue", show.colnames=FALSE, cex=1.75)
    title(main="Project Summary", family="sans", cex.main=2)
			
			
	table2<-cbind(lanes, samps)
	colnames(table2)<-c("per lane", "per sample")
	rownames(table2)<-c("Reads", "Used bases", "Average target coverage", "% loci covered to 10x", "% loci covered to 20x","% loci covered to 10x")
	par(mar=c(4,4,4,4))
	textplot(table2, rmar=1, col.rownames="dark blue", cex=1.25)
	title(main="Bases Summary", family="sans", cex.main=1.75)
	   	   	
  	table3<-cbind(summary2)
  	rownames(table3)<-c("Sequencer", "Used lanes", "Used lanes per sample", "Lane pariteies", "Read legnths")
	par(mar=c(4,4,4,4))
	textplot(table3, rmar=1, col.rownames="dark blue", show.colnames=FALSE, cex=1.25)
	title(main="Sequencing Summary", family="sans", cex.main=1.75)

	   	   			
	table4<-cbind(SNPS, titvs)
	rownames(table4)<-c("All SNPs", "Known SNPs", "Novel SNPs")
	colnames(table4)<-c("SNPs Found", "Ti/Tv")
	textplot(table4, rmar=1, col.rownames="dark blue", cex=1.25)
	title(main="Variant Summary", family="sans", cex.main=1.75)
	
	}

fingerprints<-function(lanetable, sample_sets){
	attach(lanetable)
	
	#define layout
	layout(matrix(c(1,2,3), ncol=1, nrow=3, byrow=TRUE), heights=c(1, 3,2), respect=FALSE)
    
    #prep for title bar
    title=paste(sample_sets, ": Fingerprint Status", sep="")
    drop<-read.jpeg("adprdrop.jpg")
    
    #plot title bar
    par(mar=c(0,0,0,0))
    plot(drop)
    text(100, 40, title, family="serif", adj=c(0,0), cex=3, col=gray(.25))
    
    #prep for FP plot	
    badsnps<-union(which(FP_CONFIDENT_MATCHING_SNPS<15), which(FP_CONFIDENT_MATCHING_SNPS<15))
	colors<-c(rep("Blue", length(FP_CONFIDENT_CALLS)))
    colors[badsnps]<-"Red"
    ticks<-c(match(unique(Flowcell), Flowcell) )
  	ys=rep(c(0, max(SNP_TOTAL_SNPS, na.rm=TRUE)*1.04, max(SNP_TOTAL_SNPS, na.rm=TRUE)*1.04, 0, 0), ceiling(length(ticks)/2))
    shader<-ticks[c(rep(c(1,1,2,2,1), ceiling(length(ticks)/2))+sort(rep(seq(0, length(ticks),by=2), 5)))]-0.5
    if((length(ticks)%%2 > 0)){
        shader[(length(shader)-2):(length(shader)-1)]<-length(Flowcell)+0.5
       	}
    shader<-na.omit(shader)		
    
    #plot FP plot		
    par(mar=c(10, 6, 8, 3))
    plot(1:length(FP_CONFIDENT_MATCHING_SNPS), FP_CONFIDENT_MATCHING_SNPS, pch=NA, ylim=c(0,24), ylab="Fingerprint calls", xlab="", xaxt="n", col=colors, main="Fingerprint Calling and Matching Sorted by Flowcell", cex.main=2) 
	axis(side=3, at=c(1:length(Flowcell)), labels=Lane[order(Flowcell)], cex.axis=0.5, padj=1,tick=FALSE)
	axis(side=1, at=c(ticks), labels=sort(unique(Flowcell)), tick=FALSE, las=2)
	mtext("Lane",side=3, cex=.75, line=1.5)
	mtext("Flowcell",side=1, cex=1.25, line=8)
    polygon(shader, ys, border="black", lty=0, col="gray")
    points(1:length(FP_CONFIDENT_MATCHING_SNPS), FP_CONFIDENT_MATCHING_SNPS, pch=4, col=colors)
	points(1:length(FP_CONFIDENT_MATCHING_SNPS), FP_CONFIDENT_CALLS, pch=3, col=colors)
    if(length(badsnps)>0){
    	legend("bottomright", legend=c("Confident calls at fingerprint sites by lane", "Confident matching calls at fingerprint sites by lane", "Confident calls in bad lanes", "Confident matching calls in bad lanes", "All Confident calls match fingerprint sites"), pch=c(4,3,4,3,8), col=c("Blue", "Blue", "Red", "Red", "Black" ), bg="White")
       	mtext("Some problematic fingerprint sites", side=3)
  		}else{
  	    	legend("bottomright", legend=c("Confident calls at fingerprint sites by lane", "Confident matching calls at fingerprint sites by lane", "All Confident calls match fingerprint sites"), pch=c(4, 3, 8), col="Blue", bg="White")
  	      	}

	#plot some summary of FP stuff
    textplot("Some summary of Fingerprint problems will go here ", valign="top", family="sans")

	detach(lanetable)
	}

snps_called<-function(lanetable, sample_sets){  
	attach(lanetable)  		
  	
  	#define layout for this page		
	layout(matrix(c(1,1,2, 3, 4,4), ncol=2, nrow=3, byrow=TRUE), widths = c(3,1), heights=c(1, 3,2), respect=FALSE)
    
    #prep for title bar
    title=paste(sample_sets, ": SNPs Called by Lane", sep="")
    drop<-read.jpeg("adprdrop.jpg")
    
    #plot title bar
    par(mar=c(0,0,0,0))
    plot(drop)
    text(100, 40, title, family="serif", adj=c(0,0), cex=3, col=gray(.25))

    #prep for snp plot
    ticks<-c(match(unique(Flowcell), sort(Flowcell)) )
  	ys=rep(c(min(SNP_TOTAL_SNPS, na.rm=TRUE), max(SNP_TOTAL_SNPS, na.rm=TRUE)*1.04, max(SNP_TOTAL_SNPS, na.rm=TRUE)*1.04, min(SNP_TOTAL_SNPS, na.rm=TRUE), min(SNP_TOTAL_SNPS, na.rm=TRUE)), ceiling(length(ticks)/2))
	shader<-ticks[c(rep(c(1,1,2,2,1), ceiling(length(ticks)/2))+sort(rep(seq(0, length(ticks),by=2), 5)))]-0.5
    if((length(ticks)%%2 > 0)){
        shader[(length(shader)-2):(length(shader)-1)]<-length(Flowcell)+0.5
       	}
    shader<-na.omit(shader)
    cols<-rep("blue", length(SNP_TOTAL_SNPS))
    cols[which(SNP_TOTAL_SNPS %in% boxplot.stats(SNP_TOTAL_SNPS)$out)]<-"red"
    
	#plot snp plot
	par(ylog=TRUE, mar=c(10, 6, 4, 0))
    plot(1:length(SNP_TOTAL_SNPS), SNP_TOTAL_SNPS[order(Flowcell)],xlab="", 
    	ylab="SNPs Called", 
    	ylim = c(min(SNP_TOTAL_SNPS, na.rm=TRUE), max(SNP_TOTAL_SNPS, na.rm=TRUE)), 
    	xaxt="n", 
    	pch=NA)
	title(main="SNPs Called in Each Lane sorted by Flowcell", line=3, cex=1.5)
	axis(side=3, at=c(1:length(Flowcell)), labels=Lane[order(Flowcell)], cex.axis=0.5, padj=1,tick=FALSE)
	axis(side=1, at=c(ticks), labels=sort(unique(Flowcell)), tick=FALSE, las=2)
	mtext("Lane",side=3, cex=.75, line=1.5)
	mtext("Flowcell",side=1, cex=1.25, line=8)
	polygon(shader, ys, border="black", lty=0, col="gray")
    points(1:length(SNP_TOTAL_SNPS), SNP_TOTAL_SNPS, col=cols, pch=19)
    if(length(boxplot.stats(SNP_TOTAL_SNPS)$out)>0){
    	legend("topright", legend=c("Normal SNP Call Counts", "Outlier SNP Call Counts"), pch=19, col=c("Blue", "red"), bg="White")
        }
    
    #plot boxplot
    par(ylog=TRUE, mar=c(10, 0, 4, 2))
    boxplot(SNP_TOTAL_SNPS, main="SNPs Called in Lane", ylab="", yaxt="n",  ylim = c(min(SNP_TOTAL_SNPS, na.rm=TRUE), max(SNP_TOTAL_SNPS, na.rm=TRUE)), ylog=TRUE)
	if(length(boxplot.stats(SNP_TOTAL_SNPS)$out)==0){
        mtext("No outliers",  side=1, line=4)
    	}else{
	    	mtext(paste("Outlier SNP call counts in ", length(boxplot.stats(SNP_TOTAL_SNPS)$out), "lanes"), side=1, line=4)
    	}
    
    #Plot variant summary below	
    textplot("Variant Summary will go here", valign="top", family="sans")
    
    detach(lanetable)
    }

titvsamp<-function(metricsbysamp){
	attach(titv)
	
	#define layout
	layout(matrix(c(1,2,3), ncol=1, nrow=3, byrow=TRUE), heights=c(1, 3,2), respect=FALSE)
    
    #prep for title bar
    title=paste(sample_sets, ": Ti/Tv Ratio by Sample", sep="")
    drop<-read.jpeg("adprdrop.jpg")
    
    #plot title bar
    par(mar=c(0,0,0,0))
    plot(drop)
    text(100, 40, title, family="serif", adj=c(0,0), cex=3, col=gray(.25))
    
    #prep for titv graph 
    boxplot.stats(TiTvRatio[which(filter_name=="filtered")])$stats[5]->min
    shade<-which(sort(TiTvRatio[which(novelty_name=="novel" & filter_name=="called")], decreasing=TRUE)<min)-.5
    
	#plot titv graph
	par(mar=c(9, 5, 4, 2))
	plot(seq(1:length(unique(row))), sort(TiTvRatio[which(novelty_name=="novel" & filter_name=="called")], decreasing=TRUE), 
		xaxt="n", 
		main="Ti/Tv for Novel and Known SNP calls", 
		ylab="Ti/Tv", 
		xlab="", 
		col="red",
		cex.main=2,
		cex.lab=1.25,
		cex.axis=1,
		pch=1)
	polygon(c(min(shade), min(shade), max(shade)+5, max(shade)+5, min(shade)), c(par()$xaxp[1:2], par()$xaxp[2:1], par()$xaxp[1]), col="gray", lty=0)
	points(seq(1:length(unique(row))), c(TiTvRatio[which(novelty_name=="known" & filter_name=="called")])[order(TiTvRatio[which(novelty_name=="novel" & filter_name=="called")], decreasing=TRUE)], pch=1, col="blue")
	axis(side=1, at=c(1:length(unique(row))), labels=unique(row)[order(TiTvRatio[which(novelty_name=="novel" & filter_name=="called")], decreasing=TRUE)], tick=FALSE, hadj=1, las=2, cex=1.25)
	mtext("Samples Sorted by Novel Ti/Tv", side=1, cex=1., line = 6)
	abline(a=mean(TiTvRatio[which(novelty_name=="all" & filter_name=="called")]),b=0)
	if(length(shade)<1){
		legend("topright", legend=c("Known Variants", "Novel Variants", "Mean Ti/Tv for all variants"), col=c("blue", "red", "black"), pch=c(1,1,NA), lty=c(0, 0, 1), xjust=0.5, cex=1.25, adj=c(-20, 0))
		}else{
			points(shade+.5, sort(TiTvRatio[which(novelty_name=="novel" & filter_name=="called")], decreasing=TRUE)[shade], pch=4, col="red")
			legend("top", legend=c("Known Variants", "Novel Variants (normal values)", "Novel Variants (low values)","Mean Ti/Tv for all called variants"), col=c("blue", "red", "red", "black"), pch=c(1,1,4,NA), lty=c(0, 0, 0, 1), xjust=0.5, bty="n", cex=1.25, adj=c(0, 0))
			}
	
	 #Plot variant summary below	
    par(mar=c(2, 2, 2, 2))
    textplot("Lower TiTv indicates potentially higher false positive rates.\nTi/Tv ratios within the 95% confidence interval of the distribution of Ti/Tv ratios for filtered calls are indicated by gray shading.\nSomething Else will go here too", valign="top", family="sans")
    
			
	detach(titv)

	}

#functionalclasses<-function(countfunctclasses){}

errorratepercycle<-function(erpc){

	#define layout
	layout(matrix(c(1,2,3), ncol=1, nrow=3, byrow=TRUE), heights=c(1, 3,2), respect=FALSE)
    
    #prep for title bar
    title=paste(sample_sets, ": Error Rate Per Cycle", sep="")
    drop<-read.jpeg("adprdrop.jpg")
    
    #plot title bar
    par(mar=c(0,0,0,0))
    plot(drop)
    text(100, 40, title, family="serif", adj=c(0,0), cex=3, col=gray(.25))
   
	#prep for erprp graph
	crazies<-which(errpercycle[nrow(errpercycle),]>0.3) #this can be changed to any kind of filter for particular lanes
	colors<-rainbow(ncol(errpercycle), s=0.5, v=0.5)
	colors[crazies]<-rainbow(length(crazies))
	weights<-rep(1, ncol(errpercycle))
	weights[crazies]<-2
	
	#plot erprp graph
	par(mar=c(6, 6, 3, 2))
	matplot(errpercycle, 
		type="l", 
		lty="solid", 
		col=colors, 
		lwd=weights,  
		main="Error Rate per Read Position", 
		ylab="Error Rate", 
		xlab="Cycle/Read Position", 
		log="y", 
		cex.main=2,
		cex.lab=1.5,
		cex.axis=1.25,
		)
	if(length(crazies)>0){
		legend("topleft", title="Unusual Lanes", legend=colnames(errpercycle)[crazies], lty="solid", lwd=2, col=colors[crazies], xjust=0.5)
		}else{
		mtext("No unusual lanes.", 1, line=6, cex=1.25)
		}
		
	 #Plot variant summary below	
    textplot("Something related will go here", valign="top", family="sans")
	
	}

depth_target<-function(DOC){
	
		#define layout
	layout(matrix(c(1,2), ncol=1, nrow=2, byrow=TRUE), heights=c(1, 5), respect=FALSE)
    
    #prep for title bar
    title=paste(sample_sets, ": Depth of Coverage By Target", sep="")
    drop<-read.jpeg("adprdrop.jpg")
    
    #plot title bar
    par(mar=c(0,0,0,0))
    plot(drop)
    text(100, 40, title, family="serif", adj=c(0,0), cex=1.75, col=gray(.25))

	colnames(DOC)->cols
	apply(DOC[,grep("mean", cols)], 1, median)->medianofmeans
	apply(DOC[,grep("mean", cols)], 1, quantile, probs=3/4)->q3s
	apply(DOC[,grep("mean", cols)], 1, quantile, probs=1/4)->q1s

	par(ylog=FALSE, mar=c(5, 5, 4, 2))
	plot(c(1:3122),sort(medianofmeans, decreasing=TRUE), type="l",log="y",ylab="Coverage", xlab="",xaxt="n", main="Coverage Across All Targets", lwd=2, cex.main=2.5, cex.lab=1.5, cex.axis=1.25)
	mtext("Targets sorted by median avereage coverage across sample", side=1, line=1, cex=1.5)
	abline(h=10, lty="dashed", lwd=3)
	lines(c(1:3122),q3s[order(medianofmeans, decreasing=TRUE)], col="dark blue")
	lines(c(1:3122),q1s[order(medianofmeans, decreasing=TRUE)], col="dark blue")
	legend(c(0, 20), legend="10x coverage", box.lty=0, lwd=3, lty="dashed")
	legend("bottomleft", legend=c("Median average target coverage across all samples", "First and third quartiles of average target across all sample"), box.lty=0, lwd=c(1,2), col=c("black", "dark blue"), lty="solid")

		
	#define layout
	layout(matrix(c(1,2), ncol=1, nrow=2, byrow=TRUE), heights=c(1,5), respect=FALSE)
    
    #prep for title bar
    title=paste(sample_sets, ": Depth of Coverage For Poorly Covered Targets", sep="")
    drop<-read.jpeg("adprdrop.jpg")
    
    #plot title bar
    par(mar=c(0,0,0,0))
    plot(drop)
    text(100, 40, title, family="serif", adj=c(0,0), cex=1.25, col=gray(.25))		
	yuck<-DOC[which(medianofmeans<10),grep("mean", cols)]
	yuck<-yuck+0.1
	par(mar=c(17, 4, 4, 2))
	boxplot(t(yuck[order(medianofmeans[which(medianofmeans<10)], decreasing=TRUE),]),log="y", yaxt="n", xaxt="n", cex.lab=1.15, cex.axis=1.05, ylab="Average coverage accross all samples", main="Targets with low coverage accross samples")

	axis(2, at=axTicks(2)+c(0, rep(0.1, length(axTicks(2))-1)), labels=c(0.0, axTicks(2)[2:length(axTicks(2))]), cex.axis=0.75)
	mtext("Target", side=1, line=15, cex=1.5)
	axis(1, at=c(1:length(which(medianofmeans<10))), labels=rownames(DOC[which(medianofmeans<10),])[order(medianofmeans[which(medianofmeans<10)])], las=2, cex.axis=1.15)
	}

depth_sample<-function(DOC2){
	
	#define layout
 	layout(matrix(c(1,2), ncol=1, nrow=2, byrow=TRUE), heights=c(1,5), respect=FALSE)
     
    #prep for title bar
    title=paste(sample_sets, ": Mean Depth of Coverage per Base by Sample", sep="")
    drop<-read.jpeg("adprdrop.jpg")
     
    #plot title bar
    par(mar=c(0,0,0,0))
    plot(drop)
    text(100, 40, title, family="serif", adj=c(0,0), cex=1.25, col=gray(.25))	
    #prep for bysample 
    means<-c(sort(DOC2[which(DOC2[,2]<250),2]), rep(250, (length(which(DOC2[,2]>=250))-1)))
    types<-rep(20, length(means))
    cols<-rep("black", length(means))
    types[which(means==250)]<-8
    cols[which(means==250)]<-"red"
 
 	#plot doc by sample
 	
 	par(mar=c(10, 4, 4, 2))
 	plot(means, ylim=c(0, 250), xaxt="n", col=cols, pch=types, xlab="", ylab="Depth of Coverage")	
> 	axis(1, at=c(1:(nrow(DOC2)-1)), labels=c(rownames(DOC2[which(DOC2[,2]<250),])[order(DOC2[which(DOC2[,2]<250),2])], rownames(DOC2[which(DOC2[,2]>=250),])[order(which(DOC2[,2]>=250))][1:(length(which(DOC2[,2]>=250))-1)]), las=2)
> 	mtext("Samples", side=1, line=7, cex=1.25)

	
	}



datapuller<-function(setname){
	#library(yaml)

	strsplit(setname, ".")[1]->projectname


	lanes<-read.delim(paste(projectname, "_lanes.txt", sep=""), header=TRUE)
	samps<-read.delim(paste(projectname, "_samps.txt", sep=""), header=TRUE)
	#doct<-read.delim(paste(setname, "depth.sample_interval_summary", sep=""), header=TRUE, row.names=1)
	#docs<-read.delim(paste(setname, ".depth.sample_summary", sep=""), header=TRUE, row.names=1)
	#eval<-read.csv(paste(setname, "eval.CountFunctionalClasses", sep=""), skip=1)
	titv<-read.csv(paste(setname,  ".eval.SimpleMetricsBySample.csv", sep=""), skip=1)
	#erprp<-read.delim(paste(setname, ".erprp", sep=""))
	
	colnames(lanes)<-c('Initiative','Project','GSSR.ID','External.ID','WR.ID','Flowcell','Lane','Lane.Type','Library','AL_TOTAL_READS','AL_PF_READS','AL_PCT_PF_READS','AL_PF_NOISE_READS','AL_PF_READS_ALIGNED','AL_PCT_PF_READS_ALIGNED','AL_PF_HQ_ALIGNED_READS','AL_PF_HQ_ALIGNED_BASES','AL_PF_HQ_ALIGNED_Q20_BASES','AL_PF_HQ_MEDIAN_MISMATCHES','AL_MEAN_READ_LENGTH','AL_READS_ALIGNED_IN_PAIRS','AL_PCT_READS_ALIGNED_IN_PAIRS','AL_BAD_CYCLES','AL_PCT_STRAND_BALANCE','DUP_UNPAIRED_READS_EXAMINED','DUP_READ_PAIRS_EXAMINED','DUP_UNMAPPED_READS','DUP_UNPAIRED_READ_DUPLICATES','DUP_READ_PAIR_DUPLICATES','DUP_PERCENT_DUPLICATION','DUP_ESTIMATED_LIBRARY_SIZE','HS_BAIT_SET','HS_GENOME_SIZE','HS_LIBRARY_SIZE','HS_BAIT_TERRITORY','HS_TARGET_TERRITORY','HS_BAIT_DESIGN_EFFICIENCY','HS_TOTAL_READS','HS_PF_READS','HS_PF_UNIQUE_READS','HS_PCT_PF_READS','HS_PCT_PF_UQ_READS','HS_PCT_PF_UQ_READS_ALIGNED','HS_PF_UQ_READS_ALIGNED','HS_PF_UQ_BASES_ALIGNED','HS_ON_BAIT_BASES','HS_NEAR_BAIT_BASES','HS_OFF_BAIT_BASES','HS_ON_TARGET_BASES','HS_PCT_SELECTED_BASES','HS_PCT_OFF_BAIT','HS_ON_BAIT_VS_SELECTED','HS_MEAN_BAIT_COVERAGE','HS_MEAN_TARGET_COVERAGE','HS_FOLD_ENRICHMENT','HS_ZERO_CVG_TARGETS_PCT','HS_FOLD_80_BASE_PENALTY','HS_PCT_TARGET_BASES_2X','HS_PCT_TARGET_BASES_10X','HS_PCT_TARGET_BASES_20X','HS_PCT_TARGET_BASES_30X','HS_PENALTY_10X','HS_PENALTY_20X','HS_PENALTY_30X','SNP_TOTAL_SNPS','SNP_PCT_DBSNP','SNP_NUM_IN_DBSNP','Lane.IC.Matches','Lane.IC.PCT.Mean.RD1.Err.Rate','Lane.IC.PCT.Mean.RD2.Err.Rate','FP_PANEL_NAME','FP_PANEL_SNPS','FP_CONFIDENT_CALLS','FP_CONFIDENT_MATCHING_SNPS','FP_CONFIDENT_CALLED_PCT','FP_CONFIDENT_MATCHING_SNPS_PCT','LPCNCRD_REFERENCE','LPCNCRD_NON_REFERENCE','LPCNCRD_PCT_CONCORDANCE')
	
	files<-list(c(lanes, samps, doct, docs, eval, titv, erprp))
	
	return(files)
	}


runner<-function(basename, desc1, desc2){
	datapuller(basename)->tables
	attach(tables)
	
	
	
	pdf(paste(basename, ".pdf", sep=""), width=22, height=15,pointsize=24)

	tearsheet(lanes, samps, titv, desc1, desc1)	
	fingerprints(lanes)
	snps_called(lanes)
	titvsamp(titv)
	#functionalclasses(eval)
	#errorratepercycle(erprp)
	#depth_target(doct)
	#depth_sample(docs)
	
	dev.off()
	detach(tables)
	}

if(length(commandArgs(TRUE))>0){
	runner(commandArgs(TRUE))
	}	
	



