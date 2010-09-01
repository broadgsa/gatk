#Before executing this file, save squid files as csv, then as tab deliminated files with only the column values as the header, change the format of all cells to numbers. Assign the path to these files to "samples" and "lanes" respectively. 
#set up database stuff for firehose and picard interface
#set up so runnable by firehsoe

stuffmaker<-function(args){

lanes<-args[1]
samples<-args[2]
sample_sets<-args[3]
eval<-args[4]
titveval<-args[5]
DOCi<-args[6]
DOCs<-args[7]

library(gplots)

pdf(file=paste(sample_sets, ".pdf", sep=""), width=22, height=15, pagecentre=TRUE, pointsize=24)        


if(is.na(sample_sets)){
	print("Please specify sample set for file naming and press enter.")
	scan("stdin", what="character",n=1)->sample_sets
	print("Thanks!")
	}

	if(is.na(lanes) == FALSE && is.na(samples)==FALSE){
		#this makes a table & graphs using Picard data
		
		if(typeof(lanes)=="character"){
			read.delim(file=lanes, header= TRUE)->bylane;
                        colnames(bylane)<-c('Initiative','Project','GSSR.ID','External.ID','WR.ID','Flowcell','Lane','Lane.Type','Library','AL_TOTAL_READS','AL_PF_READS','AL_PCT_PF_READS','AL_PF_NOISE_READS','AL_PF_READS_ALIGNED','AL_PCT_PF_READS_ALIGNED','AL_PF_HQ_ALIGNED_READS','AL_PF_HQ_ALIGNED_BASES','AL_PF_HQ_ALIGNED_Q20_BASES','AL_PF_HQ_MEDIAN_MISMATCHES','AL_MEAN_READ_LENGTH','AL_READS_ALIGNED_IN_PAIRS','AL_PCT_READS_ALIGNED_IN_PAIRS','AL_BAD_CYCLES','AL_PCT_STRAND_BALANCE','DUP_UNPAIRED_READS_EXAMINED','DUP_READ_PAIRS_EXAMINED','DUP_UNMAPPED_READS','DUP_UNPAIRED_READ_DUPLICATES','DUP_READ_PAIR_DUPLICATES','DUP_PERCENT_DUPLICATION','DUP_ESTIMATED_LIBRARY_SIZE','HS_BAIT_SET','HS_GENOME_SIZE','HS_LIBRARY_SIZE','HS_BAIT_TERRITORY','HS_TARGET_TERRITORY','HS_BAIT_DESIGN_EFFICIENCY','HS_TOTAL_READS','HS_PF_READS','HS_PF_UNIQUE_READS','HS_PCT_PF_READS','HS_PCT_PF_UQ_READS','HS_PCT_PF_UQ_READS_ALIGNED','HS_PF_UQ_READS_ALIGNED','HS_PF_UQ_BASES_ALIGNED','HS_ON_BAIT_BASES','HS_NEAR_BAIT_BASES','HS_OFF_BAIT_BASES','HS_ON_TARGET_BASES','HS_PCT_SELECTED_BASES','HS_PCT_OFF_BAIT','HS_ON_BAIT_VS_SELECTED','HS_MEAN_BAIT_COVERAGE','HS_MEAN_TARGET_COVERAGE','HS_FOLD_ENRICHMENT','HS_ZERO_CVG_TARGETS_PCT','HS_FOLD_80_BASE_PENALTY','HS_PCT_TARGET_BASES_2X','HS_PCT_TARGET_BASES_10X','HS_PCT_TARGET_BASES_20X','HS_PCT_TARGET_BASES_30X','HS_PENALTY_10X','HS_PENALTY_20X','HS_PENALTY_30X','SNP_TOTAL_SNPS','SNP_PCT_DBSNP','SNP_NUM_IN_DBSNP','Lane.IC.Matches','Lane.IC.PCT.Mean.RD1.Err.Rate','Lane.IC.PCT.Mean.RD2.Err.Rate','FP_PANEL_NAME','FP_PANEL_SNPS','FP_CONFIDENT_CALLS','FP_CONFIDENT_MATCHING_SNPS','FP_CONFIDENT_CALLED_PCT','FP_CONFIDENT_MATCHING_SNPS_PCT','LPCNCRD_REFERENCE','LPCNCRD_NON_REFERENCE','LPCNCRD_PCT_CONCORDANCE')
                        }else{
		 		lanes->bylane
                                colnames(bylane)<-c('Initiative','Project','GSSR.ID','External.ID','WR.ID','Flowcell','Lane','Lane.Type','Library','AL_TOTAL_READS','AL_PF_READS','AL_PCT_PF_READS','AL_PF_NOISE_READS','AL_PF_READS_ALIGNED','AL_PCT_PF_READS_ALIGNED','AL_PF_HQ_ALIGNED_READS','AL_PF_HQ_ALIGNED_BASES','AL_PF_HQ_ALIGNED_Q20_BASES','AL_PF_HQ_MEDIAN_MISMATCHES','AL_MEAN_READ_LENGTH','AL_READS_ALIGNED_IN_PAIRS','AL_PCT_READS_ALIGNED_IN_PAIRS','AL_BAD_CYCLES','AL_PCT_STRAND_BALANCE','DUP_UNPAIRED_READS_EXAMINED','DUP_READ_PAIRS_EXAMINED','DUP_UNMAPPED_READS','DUP_UNPAIRED_READ_DUPLICATES','DUP_READ_PAIR_DUPLICATES','DUP_PERCENT_DUPLICATION','DUP_ESTIMATED_LIBRARY_SIZE','HS_BAIT_SET','HS_GENOME_SIZE','HS_LIBRARY_SIZE','HS_BAIT_TERRITORY','HS_TARGET_TERRITORY','HS_BAIT_DESIGN_EFFICIENCY','HS_TOTAL_READS','HS_PF_READS','HS_PF_UNIQUE_READS','HS_PCT_PF_READS','HS_PCT_PF_UQ_READS','HS_PCT_PF_UQ_READS_ALIGNED','HS_PF_UQ_READS_ALIGNED','HS_PF_UQ_BASES_ALIGNED','HS_ON_BAIT_BASES','HS_NEAR_BAIT_BASES','HS_OFF_BAIT_BASES','HS_ON_TARGET_BASES','HS_PCT_SELECTED_BASES','HS_PCT_OFF_BAIT','HS_ON_BAIT_VS_SELECTED','HS_MEAN_BAIT_COVERAGE','HS_MEAN_TARGET_COVERAGE','HS_FOLD_ENRICHMENT','HS_ZERO_CVG_TARGETS_PCT','HS_FOLD_80_BASE_PENALTY','HS_PCT_TARGET_BASES_2X','HS_PCT_TARGET_BASES_10X','HS_PCT_TARGET_BASES_20X','HS_PCT_TARGET_BASES_30X','HS_PENALTY_10X','HS_PENALTY_20X','HS_PENALTY_30X','SNP_TOTAL_SNPS','SNP_PCT_DBSNP','SNP_NUM_IN_DBSNP','Lane.IC.Matches','Lane.IC.PCT.Mean.RD1.Err.Rate','Lane.IC.PCT.Mean.RD2.Err.Rate','FP_PANEL_NAME','FP_PANEL_SNPS','FP_CONFIDENT_CALLS','FP_CONFIDENT_MATCHING_SNPS','FP_CONFIDENT_CALLED_PCT','FP_CONFIDENT_MATCHING_SNPS_PCT','LPCNCRD_REFERENCE','LPCNCRD_NON_REFERENCE','LPCNCRD_PCT_CONCORDANCE')
		 		}
		 if(typeof(samples)=="character"){
    		read.delim(file=samples, header= TRUE)->bysample; 
			}else{
		 		samples->bysample
		 		}
		 		
		#Calc by lane metrics
		sdlane<-rep("NA", 6)
		meanlane<-sdlane
	   attach(bylane);
	   	
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
       
   		 #makes a plot of the number of SNPS called per lane
			
  
  		ticks<-c(match(unique(Flowcell), sort(Flowcell)) )
  		ys=rep(c(min(SNP_TOTAL_SNPS, na.rm=TRUE)*0.96, max(SNP_TOTAL_SNPS, na.rm=TRUE)*1.04, max(SNP_TOTAL_SNPS, na.rm=TRUE)*1.04, min(SNP_TOTAL_SNPS, na.rm=TRUE)*0.96, min(SNP_TOTAL_SNPS, na.rm=TRUE)*0.96), ceiling(length(ticks)/2))
  		
  		defaults<-par(no.readonly = TRUE)
  			
		layout(matrix(c(1,1 , 2), 1, 3, byrow=FALSE), respect=TRUE)
    	par(mar=c(10, 6, 4, 8))
    	plot(1:length(SNP_TOTAL_SNPS), SNP_TOTAL_SNPS[order(Flowcell)],xlab="", ylab="SNPs Called in Lane", ylim = c(min(SNP_TOTAL_SNPS, na.rm=TRUE), max(SNP_TOTAL_SNPS, na.rm=TRUE)), xaxt="n", pch=NA)
		title(main=paste(sample_sets, ": SNPs Called in Each Lane sorted by Flowcell", sep=""), line=3, cex=1.25)
		axis(side=3, at=c(1:length(Flowcell)), labels=Lane[order(Flowcell)], cex.axis=0.5, padj=1,tick=FALSE)
		axis(side=1, at=c(ticks), labels=sort(unique(Flowcell)), tick=FALSE, las=2)
		mtext("Lane",side=3, cex=.75, line=1.5)
				mtext("Flowcell",cex=.75,side=1, line=8)

        shader<-ticks[c(rep(c(1,1,2,2,1), ceiling(length(ticks)/2))+sort(rep(seq(0, length(ticks),by=2), 5)))]-0.5
        if((length(ticks)%%2 > 0)){
        	shader[(length(shader)-2):(length(shader)-1)]<-length(Flowcell)+0.5
        	}
        shader<-na.omit(shader)
        polygon(shader, ys, border="black", lty=0, col="gray")
        cols<-rep("blue", length(SNP_TOTAL_SNPS))
        cols[which(SNP_TOTAL_SNPS %in% boxplot.stats(SNP_TOTAL_SNPS)$out)]<-"red"
        points(1:length(SNP_TOTAL_SNPS), SNP_TOTAL_SNPS, col=cols, pch=19)
        if(length(boxplot.stats(SNP_TOTAL_SNPS)$out)>0){
          legend("topright", legend=c("Normal SNP Call Counts", "Outlier SNP Call Counts"), pch=19, col=c("Blue", "red"), bg="White")
        }
           
     	boxplot(SNP_TOTAL_SNPS, main="SNPs Called in Lane", ylab="SNPs Called"  )


		if(length(boxplot.stats(SNP_TOTAL_SNPS)$out)==0){
        	mtext("No outliers",  side=1, line=4)
    		}else{
	            mtext(paste("Outlier SNP call counts in ", length(boxplot.stats(SNP_TOTAL_SNPS)$out), "lanes"), side=1, line=4)
    			}
        
        
    		

        
    		#makes a plot of fingerprint calls and labels them good or bad
        	par(defaults)
        	
        	
        	badsnps<-union(which(FP_CONFIDENT_MATCHING_SNPS<15), which(FP_CONFIDENT_MATCHING_SNPS<15))

		    colors<-c(rep("Blue", length(FP_CONFIDENT_CALLS)))
    		colors[badsnps]<-"Red"
    		ticks<-c(match(unique(Flowcell), Flowcell) )
  			ys=rep(c(0, 24*1.04, 24*1.04, 0, 0), ceiling(length(ticks)/2))
    		#pdf(file=paste(sample_sets, "_Fingerprints.pdf", sep=""), width=.2*length(FP_CONFIDENT_CALLS), height=.1*length(FP_CONFIDENT_CALLS))
    		par(mar=c(10, 6, 8, 3))
    		plot(1:length(FP_CONFIDENT_MATCHING_SNPS), FP_CONFIDENT_MATCHING_SNPS, pch=NA, ylim=c(0,24), ylab="Fingerprint calls", xlab="", xaxt="n", col=colors, main="Fingerprint Calling and Matching Sorted by lane") 
    		axis(side=1, at=(ticks+1), labels=unique(Flowcell), tick=FALSE, hadj=1, las=2)
        shader<-ticks[c(rep(c(1,1,2,2,1), ceiling(length(ticks)/2))+sort(rep(seq(0, length(ticks),by=2), 5)))]-0.5
        shader<-na.omit(shader)
        if((length(ticks)%%2 > 0)){
        	shader[(length(shader)-2):(length(shader)-1)]<-length(Flowcell)+0.5
        	}
        
        polygon(shader, ys, border="black", lty=0, col="gray")
    		    		points(1:length(FP_CONFIDENT_MATCHING_SNPS), FP_CONFIDENT_MATCHING_SNPS, pch=4, col=colors)

    		points(1:length(FP_CONFIDENT_MATCHING_SNPS), FP_CONFIDENT_CALLS, pch=3, col=colors)
    		
  


		    if(length(badsnps)>0){
        		legend("bottomright", legend=c("Confident calls at fingerprint sites by lane", "Confident matching calls at fingerprint sites by lane", "Confident calls in bad lanes", "Confident matching calls in bad lanes", "All Confident calls match fingerprint sites"), pch=c(4,3,4,3,8), col=c("Blue", "Blue", "Red", "Red", "Black" ), bg="White")
       			 mtext("Some problematic fingerprint sites", side=3)
  				  }else{
  	      			legend("bottomright", legend=c("Confident calls at fingerprint sites by lane", "Confident matching calls at fingerprint sites by lane", "All Confident calls match fingerprint sites"), pch=c(4, 3, 8), col=c("Blue", "Blue", "Black"), bg="White")
  	      			}

 
	   		detach(bylane) 
	   		
				         
	   	}else{
	   		print("Lane and Sample metrics file paths not provided")
	   		}	 
			meansamp<-rep("NA", 6)
			sdsamp<-meansamp

	#Calc by sample metrics   	
		   	attach(bysample);
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

			#print all of this stuff out in R. 
			summary<-c(paste(callable.target, "bases"), paste(mean.lanes.samp, "+/-", sd.lanes.samp), paste(singlelanes, "single lanes,", pairedlanes, "paired lanes"), paste(mean.mrl.samp, "+/-", sd.mrl.samp))

			samps<-paste(meansamp, c("M", "M", "x", "%", "%", "%"), " +/- ", sdsamp, c("M", "M", "x", "%", "%", "%"), sep="")

			lanes<-paste(meanlane, c("M", "M", "x", "%", "%", "%"), " +/- ", sdlane, c("M", "M", "x", "%", "%", "%"), sep="")
			
			layout(matrix(c(1,2), ncol=1), heights=c(2,3))

			table1<-cbind(summary)
			rownames(table1)<-c("Callable Target", "Used Lanes per Sample", "Parities", "Read Length")
			textplot(table1, col.rownames="blue", show.colnames=FALSE, cex=1.75)
			title(main="Sequencing Summary", family="serif", cex.main=2)
			table2<-cbind(lanes, samps)
			colnames(table2)<-c("per lane", "per sample")
			rownames(table2)<-c("Reads", "Used bases", "Average target coverage", "% loci covered to 10x", "% loci covered to 20x","% loci covered to 10x")
			textplot(table2, rmar=1, col.rownames="blue", cex=1.25)
			title(main="Bases Summary", family="serif", cex.main=1.75)
	   	   		 
	
	        	par(defaults)

	#Makes Error Rate percycle graph   		 
    if(is.na(eval)==FALSE){
    	if(typeof(eval)=="character"){
    		read.delim(eval, header=TRUE)[2:ncol(read.delim(eval, header=TRUE))]->errpercycle
			}else{
				eval->errpercycle
				}
		
				
		#pdf(paste(sample_sets, "_errorrate_per_cycle.pdf", sep=""), width=6, height=5)
	
		crazies<-which(errpercycle[75,]>0.3) #this can be changed to any kind of filter for particular lanes
	
		colors<-rainbow(ncol(errpercycle), s=0.5, v=0.5)
		colors[crazies]<-rainbow(length(crazies))
		weights<-rep(1, ncol(errpercycle))
		weights[crazies]<-2
	
		matplot(errpercycle, type="l", lty="solid", col=colors, lwd=weights,  main="Error Rate per Cycle", ylab="Error Rate", xlab="Cycle", ylim=c(0, 0.7))
	
		if(length(crazies)>0){
			legend("topleft", title="Unusual Lanes", legend=colnames(errpercycle)[crazies], lty="solid", lwd=2, col=colors[crazies], xjust=0.5)
			}else{
			legend("topleft", legend="No unusual lanes.", bty="n")
			}
		
		
    	
    	}else{
    		print("Error Rate Per Cycle file paths not provided")
    		}
	
	#Makes TI/TV known v novel graph
	if(is.na(titveval)==FALSE){
		##TODO: need ot make sure this is nice and prettified. 
		titv<-read.csv(file=titveval, skip=1)
		attach(titv)

		#pdf(file=paste(sample_sets, "_TI-TV.pdf", sep=""), width=0.2*length(unique(sample)), height=0.175*length(unique(sample)))
		par(mar=c(11, 4, 4, 2))
		plot(seq(1:length(unique(sample))), Ti.Tv[which(novelty_name=="novel" & filter_name=="called")], xaxt="n", ylim=c(1, 4), main="Ti/Tv for Novel and Known SNP calls", ylab="Ti/Tv", xlab="", col="red", pch=1)

		points(seq(1:length(unique(sample))), Ti.Tv[which(novelty_name=="known" & filter_name=="called")], pch=1, col="blue")

		axis(side=1, at=(1:length(unique(sample))), labels=unique(sample), tick=FALSE, hadj=1, las=2)
		
		abline(a=mean(Ti.Tv[which(novelty_name=="all" & filter_name=="called")]),b=0)
		
		legend("bottomright", legend=c("Known Variants", "Novel 		Variants", "Mean Ti/Tv for all variants"), col=c("blue", "red", "black"), pch=c(1,1,NA_integer_), lty=c(0, 0, 1), xjust=0.5)
		mtext(line=9,"Lower Ti/Tv ratios indicate potentially increased false positive SNP rates.", side=1)
		
		
		}else{
			print("TiTV filepath not provided")
			}
			
	#Make DOC graph
	if(is.na(DOCi)==FALSE){
		#pdf(paste(sample_set, "_DOCi.pdf", sep=""), width=6, height=5)
		if(typeof(DOCi)=="character"){
			as.data.frame(read.delim(DOCi))->DOC
			}else{
			DOCi->DOCdata
			}
	
		colnames(DOC)->cols
		apply(DOC[,grep("mean", cols)], 1, median)->medianofmeans
		apply(DOC[,grep("mean", cols)], 1, quantile, probs=3/4)->q3s
		apply(DOC[,grep("mean", cols)], 1, quantile, probs=1/4)->q1s

		par(ylog=FALSE, mar=c(5, 4, 4, 2))
		plot(c(1:3122),sort(medianofmeans, decreasing=TRUE), type="l", lwd="1",log="y",ylab="Coverage", xlab="Targets sorted by median average coverage across sample",xaxt="n", main="Coverage Across All Targets")

		abline(h=10, lty="dotted")

		lines(c(1:3122),q3s[order(medianofmeans, decreasing=TRUE)])

		lines(c(1:3122),q1s[order(medianofmeans, decreasing=TRUE)])

		legend("bottomleft", "10x coverage", box.lty=0, lty="dotted")

		
		
		#pdf(paste(sample_set, "_DOCiy.pdf", sep=""), width=6, height=5)
		yuck<-DOC[which(medianofmeans<10),grep("mean", cols)]
		yuck<-yuck+0.1
		par(mar=c(16, 4, 4, 2))
		boxplot(t(yuck[order(medianofmeans[which(medianofmeans<10)], decreasing=TRUE),]),log="y", yaxt="n", xaxt="n", ylab="Average coverage accross all samples", main="Targets with low coverage accross samples")

		axis(2, at=axTicks(2)+c(0, rep(0.1, length(axTicks(2))-1)), labels=c(0.0, axTicks(2)[2:length(axTicks(2))]), cex.axis=0.75)
		mtext("Target", side=1, line=14)
		axis(1, at=c(1:length(which(medianofmeans<10))), labels=DOC[which(medianofmeans<10),1][order(medianofmeans[which(medianofmeans<10)])], las=2, cex.axis=0.75)
	
		
		
		
		}else{
			print("Depth of Coverage--intervals filepath not provided")
			}		

	if(is.na(DOCs)==FALSE){
		#pdf(paste(sample_set, "_DOCs.pdf", sep=""), width=6, height=5)
		if(typeof(DOCs)=="character"){
			as.data.frame(read.delim(DOCs))->DOC2
			}else{
			DOCs->DOCdata
			}
		par(mar=c(10, 4, 4, 2))
		boxplot(t(DOC2[,2:ncol(DOC2)]+0.1), log="y", main="Depth of Coverage by Sample", xaxt="n", yaxt="n", ylab="Coverage")

		axis(1, at=c(1:nrow(DOC2)), labels=DOC2[,1], las=2)

		axis(2, at=axTicks(2)+c(0, rep(0.1, length(axTicks(2))-1)), labels=floor(c(0.0, axTicks(2)[2:length(axTicks(2))])))

		labels=floor(c(0.0, axTicks(2)[2:length(axTicks(2))]))

		mtext("Samples", side=1, line=9)
		
		
		
		}else{
			print("Depth of Coverage--samples filepath not provided")
			}
			
		dev.off()

}
if(length(commandArgs(TRUE))>0){
	stuffmaker(commandArgs(TRUE))
	}
