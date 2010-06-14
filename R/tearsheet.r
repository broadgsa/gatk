#Before executing this file, save squid files as csv, then as tab deliminated files with only the column values as the header, change the format of all cells to numbers. Assign the path to these files to "samples" and "lanes" respectively. 

tearsheetcalc<-function(lanes, samples, sample_sets, eval1, eval2, eval3, eval4, eval5, eval6){

	read.delim(file=lanes, header= TRUE)->bylane; 
    read.delim(file=samples, header= TRUE)->bysample; 

    attach(bylane);
    callable.target<-HS_TARGET_TERRITORY[1]; 
    singlelanes<-length(which(Lane.Type=="Single"));
    pairedlanes<-length(which(Lane.Type=="Paired"));
    mean.read.lane<-signif(mean(AL_TOTAL_READS, na.rm=TRUE));
    sd.read.lane<-signif(sd(AL_TOTAL_READS, na.rm=TRUE));
    mean.ub.lane<-signif(mean(HS_ON_TARGET_BASES, na.rm=TRUE));
    sd.ub.lane<-signif(sd(HS_ON_TARGET_BASES, na.rm=TRUE));
    mean.cov.lane<-round(mean(HS_MEAN_TARGET_COVERAGE, na.rm=TRUE));
    sd.cov.lane<-round(sd(HS_MEAN_TARGET_COVERAGE, na.rm=TRUE));
    mean.10x.lane<-round(mean(HS_PCT_TARGET_BASES_10X, na.rm=TRUE));
    mean.20x.lane<-round(mean(HS_PCT_TARGET_BASES_20X, na.rm=TRUE));
    mean.30x.lane<-round(mean(HS_PCT_TARGET_BASES_30X, na.rm=TRUE));
    sd.10x.lane<-round(sd(HS_PCT_TARGET_BASES_10X, na.rm=TRUE));
    sd.20x.lane<-round(sd(HS_PCT_TARGET_BASES_20X, na.rm=TRUE));
    sd.30x.lane<-round(sd(HS_PCT_TARGET_BASES_30X, na.rm=TRUE));
            
    names<-paste(Project, " ", External.ID, "-", Lane, sep="")
        
    #makes a plot of the number of SNPS called per lane

	library(graphics)


    pdf(file=paste(sample_sets, "_SNPS.pdf", sep=""), width=0.2*length(SNP_TOTAL_SNPS), height=0.1*length(SNP_TOTAL_SNPS))        
  
	layout(matrix(c(1,1 , 2), 1, 3, byrow=FALSE), respect=TRUE)
    plot(1:length(SNP_TOTAL_SNPS), main="SNPs Called in Each Lane", SNP_TOTAL_SNPS, xlab="", ylab="SNPs Called in Lane", xaxt="n", pch=16, col="blue")
	axis(side=1, at=(1:length(SNP_TOTAL_SNPS)), labels=names, cex.axis=0.75, las=2)
        
     boxplot(SNP_TOTAL_SNPS, main="SNPs Called in Lane", ylab="SNPs Called")


    if(length(boxplot.stats(SNP_TOTAL_SNPS)$out)==0){
            mtext("No outliers",  side=1, line=4)
    }else{
            mtext(paste("Outlier SNP call counts in  ", length(boxplot.stats(SNP_TOTAL_SNPS)$out), "lanes"), side=1, line=4)
    }
        
        
    dev.off()
    
    #makes SNP plot in log scale
    
     pdf(file=paste(sample_sets, "_SNPS_log.pdf", sep=""), width=0.2*length(SNP_TOTAL_SNPS), height=0.1*length(SNP_TOTAL_SNPS))       
  
	layout(matrix(c(1,1 , 2), 1, 3, byrow=FALSE), respect=TRUE)
    plot(1:length(SNP_TOTAL_SNPS), log(SNP_TOTAL_SNPS), main="SNPs Called in Each Lane",  xlab="", ylab="Log(SNPs Called in Lane)", xaxt="n", pch=16, col="blue")
    par(ylog=TRUE)
	axis(side=1, at=(1:length(SNP_TOTAL_SNPS)), labels=names, cex.axis=0.75, las=2)
        
     boxplot(SNP_TOTAL_SNPS, main="SNPs Called in Lane", ylab="SNPs Called")


    if(length(boxplot.stats(SNP_TOTAL_SNPS)$out)==0){
            mtext("No outliers",  side=1, line=4)
    }else{
            mtext(paste("Outlier SNP call counts in  ", length(boxplot.stats(SNP_TOTAL_SNPS)$out), "lanes"), side=1, line=4)
    }
        
        
    dev.off()
    
    #makes a plot of snp calls ordered by lane
    
    pdf(file=paste(sample_sets, "_SNPS_lane.pdf", sep=""), width=0.2*length(SNP_TOTAL_SNPS), height=0.1*length(SNP_TOTAL_SNPS))       
  
	layout(matrix(c(1,1 , 2), 1, 3, byrow=FALSE), respect=TRUE)
    plot(1:length(SNP_TOTAL_SNPS), SNP_TOTAL_SNPS[order(Lane)], main="SNPs Called in Each Lane",  xlab="", ylab="Log(SNPs Called in Lane)", xaxt="n", pch=16, col="blue")
    par(ylog=TRUE)
	axis(side=1, at=(1:length(SNP_TOTAL_SNPS)), labels=names[order(Lane)], cex.axis=0.75, las=2)
        
     boxplot(SNP_TOTAL_SNPS, main="SNPs Called in Lane", ylab="SNPs Called")


    if(length(boxplot.stats(SNP_TOTAL_SNPS)$out)==0){
            mtext("No outliers",  side=1, line=4)
    }else{
            mtext(paste("Outlier SNP call counts in  ", length(boxplot.stats(SNP_TOTAL_SNPS)$out), "lanes"), side=1, line=4)
    }
    
    dev.off()    
        
    #makes a plot of fingerprint calls and labels them good or bad
        
        
    badsnps<-union(which(FP_CONFIDENT_MATCHING_SNPS<15), which(FP_CONFIDENT_MATCHING_SNPS<15))

    colors<-c(rep("Blue", length(FP_CONFIDENT_CALLS)))
    colors[badsnps]<-"Red"
        
    pdf(file=paste(sample_sets, "_Fingerprints.pdf", sep=""), width=.2*length(FP_CONFIDENT_CALLS), height=.1*length(FP_CONFIDENT_CALLS))
    par(mar=c(6, 4, 5, 4))
    plot(1:length(FP_CONFIDENT_MATCHING_SNPS), FP_CONFIDENT_MATCHING_SNPS, pch=16, ylim=c(0,24), ylab="Fingerprint calls", xlab="", xaxt="n", col=colors, main="Fingerprint Calling and Matching") 
    points(1:length(FP_CONFIDENT_MATCHING_SNPS), FP_CONFIDENT_CALLS, col=colors)
    axis(side=1, at=(1:length(FP_CONFIDENT_CALLS)), labels=names, cex.axis=0.75, las=2)

    if(length(badsnps)>0){
        legend("bottomright", legend=c("Confident calls at fingerprint sites by lane", "Confident matching calls at fingerprint sites by lane", "Confident calls in bad lanes", "Confident matching calls in bad lanes"), pch=c(1, 16, 1, 16), col=c("Blue", "Blue", "Red", "Red"))
        mtext("Some problematic fingerprint sites", side=3)
    }else{
        legend("bottomright", legend=c("Confident calls at fingerprint sites by lane", "Confident matching calls at fingerprint sites by lane"), pch=c(1, 16), col="Blue")}

    dev.off()
        
    detach(bylane);

	attach(bysample);

	mean.lanes.samp<-signif(mean(X..Lanes.included.in.aggregation, na.rm = TRUE));
	sd.lanes.samp<-signif(sd(X..Lanes.included.in.aggregation, na.rm=TRUE));
	mean.mrl.samp<-signif(mean(Mean.Read.Length, na.rm=TRUE));
	sd.mrl.samp<-signif(sd(Mean.Read.Length, na.rm=TRUE));
	mean.read.samp<-signif(mean(Total.Reads, na.rm=TRUE));
	sd.read.samp<-signif(sd(Total.Reads, na.rm=TRUE));
	mean.ub.samp<-signif(mean(On.Target.Bases..HS., na.rm=TRUE));
	sd.ub.samp<-signif(sd(On.Target.Bases..HS., na.rm=TRUE));
	mean.cov.samp<-round(mean(Mean.Target.Coverage..HS., na.rm=TRUE));
	sd.cov.samp<-round(sd(Mean.Target.Coverage..HS., na.rm=TRUE));
	mean.10x.samp<-round(mean(PCT.Target.Bases.10x..HS., na.rm=TRUE));
	mean.20x.samp<-round(mean(PCT.Target.Bases.20x..HS., na.rm=TRUE));
	mean.30x.samp<-round(mean(PCT.Target.Bases.30x..HS., na.rm=TRUE));
	sd.10x.samp<-round(sd(PCT.Target.Bases.10x..HS., na.rm=TRUE));
	sd.20x.samp<-round(sd(PCT.Target.Bases.20x..HS., na.rm=TRUE));
	sd.30x.samp<-round(sd(PCT.Target.Bases.30x..HS., na.rm=TRUE));

	detach(bysample);

	#print all of this stuff out in R. 

	print(paste("Callable Target: ", callable.target, " bases", sep=""), quote = FALSE);
	print(paste("Used Lanes per Sample: ", mean.lanes.samp, " +/- ", sd.lanes.samp, sep=""), quote=FALSE);
	print(paste("Parities: ", singlelanes, " single lanes, ", pairedlanes, " paired lanes", sep=""), quote=FALSE);
	print(paste("Read Legnths: ", mean.mrl.samp, " +/- ", sd.mrl.samp, sep=""), quote = FALSE);
	print(paste("Reads per lane: ", mean.read.lane, " +/- ", sd.read.lane, sep=""), quote = FALSE);
	print(paste("Reads per sample: ", mean.read.samp, " +/- ", sd.read.samp, sep=""), quote = FALSE);
	print(paste("Used bases per lane: ", mean.read.lane, " +/- ", sd.read.lane, sep=""), quote = FALSE);
	print(paste("Used bases per sample: ", mean.read.samp, " +/- ", sd.read.samp, sep=""), quote = FALSE)
	print(paste("Average target coverage per lane: ", mean.cov.lane, " +/- ", sd.cov.lane, sep=""), quote = FALSE);
	print(paste("Average target coverage per sample: ", mean.cov.samp, " +/- ", sd.cov.samp, sep=""), quote = FALSE);
	print(paste("% loci covered to 10x per lane: ", mean.10x.lane, "% +/- ", sd.10x.lane, "%", sep=""), quote = FALSE)
	print(paste("% loci covered to 10x per sample: ", mean.10x.samp, " +/- ", sd.10x.samp, "%", sep=""), quote = FALSE)
	print(paste("% loci covered to 20x per lane: ", mean.20x.lane, "% +/- ", sd.20x.lane, "%", sep=""), quote = FALSE)
	print(paste("% loci covered to 20x per sample: ", mean.20x.samp, "% +/- ", sd.20x.samp, "%", sep=""), quote = FALSE)
	print(paste("% loci covered to 30x per lane: ", mean.30x.lane, "% +/- ", sd.30x.lane, "%", sep=""), quote = FALSE)
	print(paste("% loci covered to 30x per sample: ", mean.30x.samp, "% +/- ", sd.30x.samp, "%", sep=""), quote = FALSE)

	#still need to figure out how to get various information from eval file into R, but whatever

	read.csv(eval1, header=TRUE)->errpercycle
	
	
	errpercycle<-matrix(c(rep("tester", 24*90 ), rep(1:8, 3*90), rep(1:90,24), runif(90*24, min=0, max=0.7)), nrow=90*24, ncol=4); colnames(errpercycle)<-c("sample_set", "plate", "lane", "cycle", "errorrate") #delete this after testing
	as.data.frame(errpercycle)->errpercycle
	attach(errpercycle)
	
	pdf(paste(sample_set, "error_rate_per_cycle.pdf", sep=""))
	names<-paste(sample)


}

lane<-"mennonite_by_lane.txt"
#insert file path here before running
sample<-"mennonite_by_sample.txt"
#insertfilepath here before running
sample_set<-"Ashuldin_mennonite"


tearsheetcalc(lane,sample,sample_set)

