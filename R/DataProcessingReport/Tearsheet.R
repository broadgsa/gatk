#New tearsheet generator
.libPaths('/humgen/gsa-pipeline/.repository/R/') 

suppressMessages(library(gplots));
suppressMessages(library(ReadImages));
suppressMessages(library(gsalib));

tearsheet<-function(){
  
def.par <- par(no.readonly = TRUE)

	#define layout
	postable<-matrix(c(1, 1,  1, 1, rep(c(2, 2, 4, 4), 5), rep(c(3, 3, 4, 4), 3), rep(c(3,3,5,5), 5), 6,7,8,9), nrow=15, ncol=4, byrow=TRUE)
	layout(postable, heights=c(1, rep(.18, 13), 2), respect=FALSE)
    
    #prep for title bar
    drop<-read.jpeg(system.file("data", "tearsheetdrop.jpg", package="gsalib"))
    
    
    #plot title bar
    par(mar=c(0,0,0,0))
    plot(drop)
    text(155, 50, cmdargs$title, family="serif", adj=c(0,0), cex=3, col=gray(.25))
      print("Title created...")


	# Project summary
	projects = paste(squids, collapse=", ");

	used_samples = nrow(settable);

	unused_samples = 0;

	sequencing_protocol = samp$Initiative[1]

	bait_design = samp$"Bait Set"[1]

	callable_target = samp$"Target Territory"[1]

	table1<-rbind(paste(used_samples," used samples/", unused_samples + used_samples," total samples", sep=""), sequencing_protocol, bait_design, callable_target)
	rownames(table1)<-c("Samples","Sequencing Initiative", "Bait Design","Callable Target")
	par(mar=c(0,0,1,0))
	textplot(table1, col.rownames="darkblue", show.colnames=FALSE, cex=1.25, valign="top")
  title(main=sprintf("Project Summary (%s)\n", projects), family="sans", cex.main=1.25, line=-1)	
  print("Project summary created...")
  # Bases summary 

	reads_per_lane_mean = format(mean(lane$"PF Reads (HS)", na.rm=TRUE), 8, 3,1, scientific=TRUE);
	reads_per_lane_sd = format(sd(lane$"PF Reads (HS)", na.rm=TRUE), 8, 3,1, scientific=TRUE);
	lanessum<-sprintf("%s +/- %s\n", reads_per_lane_mean, reads_per_lane_sd)

	used_bases_per_lane_mean = format(mean(lane$"PF HQ Aligned Q20 Bases", na.rm=TRUE),8, 3,1, scientific=TRUE);
	used_bases_per_lane_sd = format(sd(lane$"PF HQ Aligned Q20 Bases", na.rm=TRUE), 8, 3,1, scientific=TRUE);
	lanessum<-c(lanessum, sprintf("%s +/- %s\n", used_bases_per_lane_mean, used_bases_per_lane_sd));

	target_coverage_mean = mean(na.omit(lane$"Mean Target Coverage"));
	target_coverage_sd = sd(na.omit(lane$"Mean Target Coverage"));
	lanessum<-c(lanessum, sprintf("%0.2fx +/- %0.2fx\n", target_coverage_mean, target_coverage_sd));

	pct_loci_gt_10x_mean = mean(na.omit(lane$"Target Bases 10x %"));
	pct_loci_gt_10x_sd = sd(na.omit(lane$"Target Bases 10x %"));
	lanessum<-c(lanessum, sprintf("%0.2f%% +/- %0.2f%%\n", pct_loci_gt_10x_mean, pct_loci_gt_10x_sd));

	pct_loci_gt_20x_mean = mean(na.omit(lane$"Target Bases 20x %"));
	pct_loci_gt_20x_sd = sd(na.omit(lane$"Target Bases 20x %"));
	lanessum<-c(lanessum,sprintf("%0.2f%% +/- %0.2f%%\n", pct_loci_gt_20x_mean, pct_loci_gt_20x_sd));

	pct_loci_gt_30x_mean = mean(na.omit(lane$"Target Bases 30x %"));
	pct_loci_gt_30x_sd = sd(na.omit(lane$"Target Bases 30x %"));
	lanessum<-c(lanessum,sprintf("%0.2f%% +/- %0.2f%%\n", pct_loci_gt_30x_mean, pct_loci_gt_30x_sd));


	reads_per_sample_mean = format(mean(samp$"PF Reads", na.rm=TRUE), 8, 3,1, scientific=TRUE);
	reads_per_sample_sd = format(sd(samp$"PF Reads",na.rm=TRUE), 8, 3,1, scientific=TRUE);
	sampssum<-sprintf("%s +/- %s\n", reads_per_sample_mean, reads_per_sample_sd);

	used_bases_per_sample_mean = format(mean(samp$"PF HQ Aligned Q20 Bases", na.rm=TRUE),8, 3,1, scientific=TRUE);
	used_bases_per_sample_sd = format(sd(samp$"PF HQ Aligned Q20 Bases", na.rm=TRUE), 8, 3,1, scientific=TRUE);
	sampssum<-c(sampssum, sprintf("%s +/- %s\n", used_bases_per_sample_mean, 	used_bases_per_sample_sd));

	target_coverage_mean = mean(na.omit(samp$"Mean Target Coverage"));
	target_coverage_sd = sd(na.omit(samp$"Mean Target Coverage"));
	sampssum<-c(sampssum, sprintf("%0.2fx +/- %0.2fx\n", target_coverage_mean, target_coverage_sd));

	pct_loci_gt_10x_mean = mean(na.omit(samp$"Target Bases 10x %"));
	pct_loci_gt_10x_sd = sd(na.omit(samp$"Target Bases 10x %"));
	sampssum<-c(sampssum, sprintf("%0.2f%% +/- %0.2f%%\n", pct_loci_gt_10x_mean, pct_loci_gt_10x_sd));

	pct_loci_gt_20x_mean = mean(na.omit(samp$"Target Bases 20x %"));
	pct_loci_gt_20x_sd = sd(na.omit(samp$"Target Bases 20x %"));
	sampssum<-c(sampssum, sprintf("%0.2f%% +/- %0.2f%%\n", pct_loci_gt_20x_mean, pct_loci_gt_20x_sd));

	pct_loci_gt_30x_mean = mean(na.omit(samp$"Target Bases 30x %"));
	pct_loci_gt_30x_sd = sd(na.omit(samp$"Target Bases 30x %"));
	sampssum<-c(sampssum, sprintf("%0.2f%% +/- %0.2f%%\n", pct_loci_gt_30x_mean, pct_loci_gt_30x_sd));

	table2<-cbind(lanessum, sampssum)
	used_lanes = length(unique(paste(lane$Flowcell, lane$Lane)));
  if(nrow(lane)>used_lanes){
    colnames(table2)<-c("Per barcoded readgroup", "Per sample")
  }
  else{
    colnames(table2)<-c("Per lane", "Per sample")
  }
	rownames(table2)<-c("Reads", "Used bases", "Average target coverage", "% loci covered to 10x", "% loci covered to 20x","% loci covered to 30x")
	par(mar=c(0,0,1,0))
	textplot(table2, rmar=1, col.rownames="dark blue", cex=1.25, valign="top")
	title(main="Bases Summary", family="sans", cex.main=1.25, line=0)

  print("Bases summary created...")

# Sequencing summary

	instrument <- c();
	if(length(grep("AAXX", lane$Flowcell))>0){
		instrument <- c(instrument, "Illumina GA2")
		}
	if(length(grep("ABXX", lane$Flowcell))>0){
		instrument <- c(instrument, "Illumina HiSeq")
		}

	if(length(instrument)>1){
		instrument<-paste(instrument[1], instrument[2], sep=" and ")
		}	

  used_lanes = length(unique(paste(lane$Flowcell, lane$Lane)));
	unused_lanes_by_sequencing = 0; #can we get this?
	unused_lanes_by_analysis = 0;

	lanes_per_sample_mean = mean(table(lane$"External ID"), na.rm=TRUE);
	lanes_per_sample_sd = sd(table(lane$"External ID"), na.rm=TRUE);
	lanes_per_sample_median = median(table(lane$"External ID"));
	lanes_paired = length(unique(paste(subset(lane, lane$"Lane Type" == "Paired")$Flowcell, subset(lane, lane$"Lane Type" == "Paired")$Lane)));
	lanes_widowed = length(unique(paste(subset(lane, lane$"Lane Type" == "Widowed")$Flowcell, subset(lane, lane$"Lane Type" == "Widowed")$Lane)));
	lanes_single = length(unique(paste(subset(lane, lane$"Lane Type" == "Single")$Flowcell, subset(lane, lane$"Lane Type" == "Single")$Lane)));

  read_length_mean = mean(lane$"Mean Read Length (P)");
	read_length_sd = sd(lane$"Mean Read Length (P)");
	read_length_median = median(lane$"Mean Read Length (P)");


	date = sort(as.Date(lane$"Run Date", format="%d-%b-%y"));

	start_date = format(date[1], "%B %d, %Y");
	end_date = format(date[length(date)], "%B %d, %Y");
 
if(nrow(lane)>used_lanes){
    used_lanes<-paste(used_lanes, " (multiplexed; ", nrow(lane), " total barcoded readgroups)", sep="")
	}
	table3<-rbind(paste(instrument), used_lanes, sprintf("%s rejected by sequencing, %s by analysis\n", unused_lanes_by_sequencing, unused_lanes_by_analysis), sprintf("%0.1f +/- %0.1f lanes (median=%0.1f)\n", lanes_per_sample_mean, lanes_per_sample_sd, lanes_per_sample_median), sprintf("%s paired, %s widowed, %s single\n", lanes_paired, lanes_widowed, lanes_single), sprintf("%0.1f +/- %0.1f bases (median=%0.1f)\n", read_length_mean, read_length_sd, read_length_median), sprintf("\tSequencing dates: %s to %s\n", start_date, end_date))
	
  rownames(table3)<-c("Sequencer", "Used lanes", "Unused lanes","Used lanes/sample", "Lane parities", "Read lengths", "Sequencing dates")
	par(mar=c(0,0,1,0))
	textplot(table3, rmar=1, col.rownames="dark blue", show.colnames=FALSE, cex=1.25, valign="top")
	title(main="Sequencing Summary", family="sans", cex.main=1.25, line=0)

  print("Sequencing summary created...")

# Variant summary
 
	eval.counts = basiceval$CountVariants
  if("FunctionalClass" %in% colnames(eval.counts)){
    eval.counts= subset(eval.counts, FunctionalClass == "all")
  }
  if("Sample" %in% colnames(eval.counts)){
    eval.counts= subset(eval.counts, Sample == "all")
  }
  if("Filter" %in% colnames(eval.counts)){
    eval.counts= subset(eval.counts, Filter == "called")
  }
	eval.counts.all = subset(eval.counts,  Novelty == "all")$nVariantLoci;
	eval.counts.known = subset(eval.counts,Novelty == "known")$nVariantLoci;
	eval.counts.novel = subset(eval.counts, Novelty == "novel")$nVariantLoci;

	eval.titv = basiceval$TiTvVariantEvaluator
	if("FunctionalClass" %in% colnames(eval.titv)){
    eval.titv= subset(eval.titv, FunctionalClass == "all")
  }
  if("Sample" %in% colnames(eval.titv)){
    eval.titv= subset(eval.titv, Sample == "all")
  }
  if("Filter" %in% colnames(eval.titv)){
    eval.titv= subset(eval.titv, Filter == "called")
  }
  eval.titv.all = subset(eval.titv, Novelty == "all")$tiTvRatio;
	eval.titv.known = subset(eval.titv, Novelty == "known")$tiTvRatio;
	eval.titv.novel = subset(eval.titv, Novelty == "novel")$tiTvRatio;

	table4 = matrix(c(eval.counts.all, eval.counts.known, eval.counts.novel, eval.titv.all, eval.titv.known, eval.titv.novel, "3.0 - 3.2", "3.2 - 3.4", "2.7 - 3.0"), nrow=3);

	rownames(table4) = c("All", "Known", "Novel");
	colnames(table4) = c("Found", "Ti/Tv ratio", "Expected Ti/Tv ratio");

	print("Variant summary created...")

	par(mar=c(0,0,0,0))
	textplot(table4, rmar=1, col.rownames="dark blue", cex=1.25, valign="top")
	title(main="Variant Summary", family="sans", cex.main=1.25, line=-2)

	eval.bysample = SAeval$CountVariants
	eval.bysample.all = subset(eval.bysample, Novelty == "all" & Sample != "all");
	eval.bysample.known = subset(eval.bysample, Novelty == "known"& Sample != "all");
	eval.bysample.novel = subset(eval.bysample, Novelty == "novel"& Sample != "all");

  eval.bysampleTITV = SAeval$TiTvVariantEvaluator
	eval.bysampleTITV.all = subset(eval.bysampleTITV, Novelty == "all" & Sample != "all");
	eval.bysampleTITV.known = subset(eval.bysampleTITV, Novelty == "known"& Sample != "all");
	eval.bysampleTITV.novel = subset(eval.bysampleTITV, Novelty == "novel"& Sample != "all");


  eval.ac = basiceval$SimpleMetricsByAC.metrics
  if("FunctionalClass" %in% colnames(eval.titv)){
    eval.ac= subset(eval.ac, FunctionalClass == "all")
  }
  if("Sample" %in% colnames(eval.titv)){
    eval.ac= subset(eval.ac, Sample == "all")
  }
  if("Filter" %in% colnames(eval.titv)){
    eval.ac= subset(eval.ac, Filter == "called")
  }
  
  eval.ac.all = subset(eval.ac, Novelty == "all");
	eval.ac.known = subset(eval.ac, Novelty == "known");
  eval.ac.novel = subset(eval.ac, Novelty == "novel");
  
	eval.func = FCeval$CountVariants

  par(mar=c(5, 5, 4, 2) + 0.1)

  
	boxplot(eval.bysampleTITV.all$tiTvRatio, eval.bysampleTITV.known$tiTvRatio, eval.bysampleTITV.novel$tiTvRatio, main="Ti/Tv by Sample", col=c("dark gray", "blue", "red"), names=c("All", "Known", "Novel"), ylab="Ti/Tv per sample", main="",cex=1.3, cex.lab=1.3, cex.axis=1.3);

	par(mar=c(7, 5, 4, 2) + 0.1)
	ind = order(eval.bysample.all$nVariantLoci);
	plot(eval.bysample.all$nVariantLoci[ind], xlab="",pch=16, col="black", xaxt="n", cex=1.1, cex.lab=1.1, cex.axis=1.1, main="Variants per Sample", ylab="Number of variants\n(axis in log space)", bty="n", log="y",ylim=c(1, max(eval.bysample.all$nVariantLoci)));
	points(eval.bysample.known$nVariantLoci[ind], pch=16, col="blue", cex=1.3);
	points(eval.bysample.novel$nVariantLoci[ind], pch=16,col="red", cex=1.3);
	legend("bottomleft", max(eval.bysample.all$nVariantLoci)/2, c("All", "Known", "Novel"), , col=c("black", "blue", "red"), pt.cex=1.3, pch=16);
  if(nrow(samp)<25){
    axis(1, at=c(1:length(eval.bysample.all$Sample[ind])), lab=eval.bysample.all$Sample[ind], cex=.7, las=2 )
  }else{
    axis(1, at=c(1:nrow(samp)), lab=rep("", nrow(samp)), cex=0.1, las=2, lwd.ticks=0)
    title(xlab="Sample\n(too many individuals to label)")
  }

	par(mar=c(6, 5, 4, 2) + 0.1)
	plot(sort(eval.ac.all$AC), eval.ac.all$n[order(eval.ac.all$AC)], ylim=c(1, max(eval.ac$n)), col="black", type="l", lwd=2, cex=1.1, cex.lab=1.1, cex.axis=1.1, xlab="Allele count\n(axis in log space)", ylab="Number of variants\n(axis in log space)", main="Variants by Allele Count", log="xy", bty="n");
	points(sort(eval.ac.known$AC), eval.ac.known$n[order(eval.ac.known$AC)], col="blue", type="l", lwd=2);
	points(sort(eval.ac.novel$AC), eval.ac.novel$n[order(eval.ac.novel$AC)], col="red", type="l", lwd=2);
	if(nrow(samp)<25){
    legend("bottomleft", c("All", "Known", "Novel"), col=c("black", "blue", "red"), lwd=2);
	}else{
    legend("topright", c("All", "Known", "Novel"), col=c("black", "blue", "red"), lwd=2);
  }
  par(mar=c(5, 5, 4, 2) + 0.1)

  barplot(eval.func$nVariantLoci[4:nrow(eval.func)], col=c("dark gray", "blue", "red"), space=c(.2,0,0), log="y", main="Variants by Functional Class", xlab="Functional Class", ylab="Number of variants\n(axis in log space)")
  axis(1, at=c(1.5,5,8.5), lab=c("Missense", "Nonsense", "Silent"), cex=.5, tick=FALSE)
  legend("top", c("All", "Known", "Novel"), fill=c("dark gray", "blue", "red"), cex=.7);

  print("Graphs created...")
	
  print("All done!")
  par(def.par)#- reset to default
	}

