##put titles/rownames left
##make titles blue
##decrease margins below titles
## put row names in black
##put background rows in. 
##change layouts so that it looks better
##get sample numbers in correctly

.libPaths('/humgen/gsa-firehose2/pipeline/repositories/StingProduction/R/') 

suppressMessages(library(gplots));
suppressMessages(library(ReadImages));

suppressMessages(library(gsalib));
suppressMessages(library(ROracle));

cmdargs = gsa.getargs(
    list(
        yaml = list(value=NA, doc="pipeline YAML file"),
        bamlist = list(value=NA, doc="list of BAM files"),
        evalroot = list(value=NA, doc="VariantEval file"),
        tearout = list(value=NA, doc="Output path for tearsheet PDF")#,
        plotout = list(value=NA, doc="Output path for PDF")
    ),
    doc="Creates a tearsheet"
);


bamlist = scan(cmdargs$bamlist, "character");
squids <- system(paste("grep SQUID ", cmdargs$yaml, ' |grep "C..." -o',  sep=""), intern=TRUE)
indexed = c();
nonindexed = c();
for (bam in bamlist) {
    bamheader = system(paste("samtools view -H", bam), intern=TRUE);
    

    if (length(bamheader) > 0) {
        rgs = bamheader[grep("^@RG", bamheader)];

        for (rg in rgs) {
            id = grep("PU:", unlist(strsplit(rg, "\t")), value=TRUE);
            id = sub("PU:", "", id);
            id = gsub("XX......", "XX", id)
	    if (length(unlist(strsplit(id, "\\.")))==3){
	       indexed<-c(indexed, id)
	       }
	       else{
		if(length(unlist(strsplit(id, "\\.")))==2){
		   nonindexed<-c(nonindexed, id)
		   }
		   else{
		      print(id + " is a strange PU and will result in odd searches")
		      }
		}
	    }
    } else {
        print(sprintf("Could not load '%s'\n", bam));
    }
}

drv = dbDriver("Oracle");
con = dbConnect(drv, "REPORTING/REPORTING@ora01:1521/SEQPROD");

rs  = dbSendQuery(con, statement = paste("SELECT * FROM ILLUMINA_PICARD_METRICS"));
d = fetch(rs, n=-1);
dbHasCompleted(rs);
dbClearResult(rs);

rs2 = dbSendQuery(con, statement = paste("SELECT * FROM ILLUMINA_SAMPLE_STATUS_AGG"));
d2 = fetch(rs2, n=-1);
dbHasCompleted(rs2);
dbClearResult(rs2);

oraCloseDriver(drv);

squid_fclanes = sprintf("%s.%s", d$"Flowcell", d$"Lane");
squid_fclanes_indexed = sprintf("%s.%s.%s", d$"Flowcell", d$"Lane", d$"Barcode");


dproj = d[which(squid_fclanes %in% nonindexed),];
dproj = rbind(dproj, d[which(squid_fclanes_indexed %in% indexed),])

dproj = dproj[which(dproj$"Project" %in% unique(squids)),]

d2proj = d2[which(d2$"Project" %in% unique(dproj$Project) & d2$"Sample" %in% dproj$"External ID"),];



tearsheet<-function(){
	tearsheetdrop <- "~Documents/Sting/R/gsalib/data/tearsheetdrop.jpg" #put the path to the tearsheet backdrop here

	pdf(file= cmdargs$tearout, width=22, height=17, pagecentre=TRUE, pointsize=24) 
  
	#define layout
	postable<-matrix(c(1, 1, 1, 1, 1, 1, rep(c(2, 2, 2, 4, 4, 4), 5), rep(c(3, 3, 3, 4, 4, 4), 3), rep(c(3,3,3,5,5,5), 5), 6,6,6,7,7,7), nrow=15, ncol=6, byrow=TRUE)
	layout(postable, heights=c(1, rep(.18, 13), 2), respect=FALSE)
    
    
    #prep for title bar
    drop<-read.jpeg(system.file(tearsheetdrop, package="gsalib"))
    
    #plot title bar
    par(mar=c(0,0,0,0))
    plot(drop)
    text(155, 50, "testing", family="serif", adj=c(0,0), cex=3, col=gray(.25))
    

	# Project summary
	projects = paste(unique(dproj$"Project"), collapse=", ");

	used_samples = length(bamlist);

	unused_samples = 0;

	sequencing_protocol = "Hybrid selection"; #can this be extracted?

	bait_design = paste(dimnames(table(dproj$"Bait Set"))[[1]][order(table(dproj$"Bait Set"), decreasing=TRUE)], collapse=", ");

	if(nchar(bait_design)>50){
		bait_design<-strsplit(bait_design, ", ")[[1]][1]	
		}

	if(nchar(bait_design)>50){
		bait_design<-strsplit(bait_design, ".Homo")[[1]][1]
		}

	callable_target = paste(na.omit(unique(dproj$"Target Territory")), collapse=", ");

	table1<-rbind(paste(used_samples," used samples/", unused_samples + used_samples," total samples", sep=""), sequencing_protocol, bait_design, callable_target)
	rownames(table1)<-c("Samples","Sequencing Protocol", "Bait Design","Callable Target")
	par(mar=c(0,0,1,0))
	textplot(table1, col.rownames="darkblue", show.colnames=FALSE, cex=1.25, valign="top")
  title(main=sprintf("Project Summary (%s)\n", projects), family="sans", cex.main=1.25, line=-1)	
  
  # Bases summary 

	reads_per_lane_mean = format(mean(dproj$"PF Reads (HS)", na.rm=TRUE), 8, 3,1, scientific=TRUE);
	reads_per_lane_sd = format(sd(dproj$"PF Reads (HS)", na.rm=TRUE), 8, 3,1, scientific=TRUE);
	lanes<-sprintf("%s +/- %s\n", reads_per_lane_mean, reads_per_lane_sd)

	used_bases_per_lane_mean = format(mean(dproj$"PF HQ Aligned Q20 Bases", na.rm=TRUE),8, 3,1, scientific=TRUE);
	used_bases_per_lane_sd = format(sd(dproj$"PF HQ Aligned Q20 Bases", na.rm=TRUE), 8, 3,1, scientific=TRUE);
	lanes<-c(lanes, sprintf("%s +/- %s\n", used_bases_per_lane_mean, used_bases_per_lane_sd));

	target_coverage_mean = mean(na.omit(dproj$"Mean Target Coverage"));
	target_coverage_sd = sd(na.omit(dproj$"Mean Target Coverage"));
	lanes<-c(lanes, sprintf("%0.2fx +/- %0.2fx\n", target_coverage_mean, target_coverage_sd));

	pct_loci_gt_10x_mean = mean(na.omit(dproj$"Target Bases 10x %"));
	pct_loci_gt_10x_sd = sd(na.omit(dproj$"Target Bases 10x %"));
	lanes<-c(lanes, sprintf("%0.2f%% +/- %0.2f%%\n", pct_loci_gt_10x_mean, pct_loci_gt_10x_sd));

	pct_loci_gt_20x_mean = mean(na.omit(dproj$"Target Bases 20x %"));
	pct_loci_gt_20x_sd = sd(na.omit(dproj$"Target Bases 20x %"));
	lanes<-c(lanes,sprintf("%0.2f%% +/- %0.2f%%\n", pct_loci_gt_20x_mean, pct_loci_gt_20x_sd));

	pct_loci_gt_30x_mean = mean(na.omit(dproj$"Target Bases 30x %"));
	pct_loci_gt_30x_sd = sd(na.omit(dproj$"Target Bases 30x %"));
	lanes<-c(lanes,sprintf("%0.2f%% +/- %0.2f%%\n", pct_loci_gt_30x_mean, pct_loci_gt_30x_sd));


	reads_per_sample_mean = format(mean(d2proj$"PF Reads", na.rm=TRUE), 8, 3,1, scientific=TRUE);
	reads_per_sample_sd = format(sd(d2proj$"PF Reads",na.rm=TRUE), 8, 3,1, scientific=TRUE);
	samps<-sprintf("%s +/- %s\n", reads_per_sample_mean, reads_per_sample_sd);

	used_bases_per_sample_mean = format(mean(d2proj$"PF HQ Aligned Q20 Bases", na.rm=TRUE),8, 3,1, scientific=TRUE);
	used_bases_per_sample_sd = format(sd(d2proj$"PF HQ Aligned Q20 Bases", na.rm=TRUE), 8, 3,1, scientific=TRUE);
	samps<-c(samps, sprintf("%s +/- %s\n", used_bases_per_sample_mean, 	used_bases_per_sample_sd));

	target_coverage_mean = mean(na.omit(d2proj$"Mean Target Coverage"));
	target_coverage_sd = sd(na.omit(d2proj$"Mean Target Coverage"));
	samps<-c(samps, sprintf("%0.2fx +/- %0.2fx\n", target_coverage_mean, target_coverage_sd));

	pct_loci_gt_10x_mean = mean(na.omit(d2proj$"Target Bases 10x %"));
	pct_loci_gt_10x_sd = sd(na.omit(d2proj$"Target Bases 10x %"));
	samps<-c(samps, sprintf("%0.2f%% +/- %0.2f%%\n", pct_loci_gt_10x_mean, pct_loci_gt_10x_sd));

	pct_loci_gt_20x_mean = mean(na.omit(d2proj$"Target Bases 20x %"));
	pct_loci_gt_20x_sd = sd(na.omit(d2proj$"Target Bases 20x %"));
	samps<-c(samps, sprintf("%0.2f%% +/- %0.2f%%\n", pct_loci_gt_20x_mean, pct_loci_gt_20x_sd));

	pct_loci_gt_30x_mean = mean(na.omit(d2proj$"Target Bases 30x %"));
	pct_loci_gt_30x_sd = sd(na.omit(d2proj$"Target Bases 30x %"));
	samps<-c(samps, sprintf("%0.2f%% +/- %0.2f%%\n", pct_loci_gt_30x_mean, pct_loci_gt_30x_sd));

	table2<-cbind(lanes, samps)
	colnames(table2)<-c("Per lane", "Per sample")

	rownames(table2)<-c("Reads", "Used bases", "Average target coverage", "% loci covered to 10x", "% loci covered to 20x","% loci covered to 30x")
	par(mar=c(0,0,1,0))
	textplot(table2, rmar=1, col.rownames="dark blue", cex=1.25, valign="top")
	title(main="Bases Summary", family="sans", cex.main=1.25, line=0)


# Sequencing summary

	instrument <- c();
	if(length(grep("AAXX", dproj$Flowcell))>0){
		instrument <- c(instrument, "Illumina GA2")
		}
	if(length(grep("ABXX", dproj$Flowcell))>0){
		instrument <- c(instrument, "Illumina HiSeq")
		}

	if(length(instrument)>1){
		instrument<-paste(instrument[1], instrument[2], sep=" and ")
		}	

	used_lanes = nrow(dproj);
	unused_lanes_by_sequencing = 0; #can we get this?
	unused_lanes_by_analysis = 0;


	lanes_per_sample_mean = mean(table(dproj$"External ID"), na.rm=TRUE);
	lanes_per_sample_sd = sd(table(dproj$"External ID"), na.rm=TRUE);
	lanes_per_sample_median = median(table(dproj$"External ID"));
	lanes_paired = nrow(subset(dproj, dproj$"Lane Type" == "Paired"));
	lanes_widowed = nrow(subset(dproj, dproj$"Lane Type" == "Widowed"));
	lanes_single = nrow(subset(dproj, dproj$"Lane Type" == "Single"));

	read_length_mean = mean(dproj$"Mean Read Length (P)");
	read_length_sd = sd(dproj$"Mean Read Length (P)");
	read_length_median = median(dproj$"Mean Read Length (P)");

	date = dproj$"Run Date";
#	date = sub("JAN", "01", date);
#	date = sub("FEB", "02", date);
#	date = sub("MAR", "03", date);
#	date = sub("APR", "04", date);
#	date = sub("MAY", "05", date);
#	date = sub("JUN", "06", date);
#	date = sub("JUL", "07", date);
#	date = sub("AUG", "08", date);
#	date = sub("SEP", "09", date);
#	date = sub("OCT", "10", date);
#	date = sub("NOV", "11", date);
#	date = sub("DEC", "12", date);
	date = date[order(as.Date(date, format="%d-%m-%Y"))];

	start_date = date[1];
	end_date = date[length(date)];
	
	
	table3<-rbind(paste(instrument), used_lanes, sprintf("%s rejected by sequencing, %s by analysis\n", unused_lanes_by_sequencing, unused_lanes_by_analysis), sprintf("%0.1f +/- %0.1f lanes (median=%0.1f)\n", lanes_per_sample_mean, lanes_per_sample_sd, lanes_per_sample_median), sprintf("%s paired, %s widowed, %s single\n", lanes_paired, lanes_widowed, lanes_single), sprintf("%0.1f +/- %0.1f bases (median=%0.1f)\n", read_length_mean, read_length_sd, read_length_median), sprintf("\tSequencing dates: %s to %s\n", start_date, end_date))
	

  	rownames(table3)<-c("Sequencer", "Used lanes", "Unused lanes","Used lanes/sample", "Lane parities", "Read lengths", "Sequencing dates")
	par(mar=c(0,0,1,0))
	textplot(table3, rmar=1, col.rownames="dark blue", show.colnames=FALSE, cex=1.25, valign="top")
	title(main="Sequencing Summary", family="sans", cex.main=1.25, line=0)

eval = gsa.read.gatkreport(cmdargs$evalroot)


# Variant summary
##TODO: Fix this csv reader
	eval.counts = eval$CountVariants
	eval.counts.all = subset(eval.counts, Novelty == "all")$nVariantLoci;
	eval.counts.known = subset(eval.counts, Novelty == "known")$nVariantLoci;
	eval.counts.novel = subset(eval.counts, Novelty == "novel")$nVariantLoci;

	eval.titv = eval$TiTvVariantEvaluator
	eval.titv.all = subset(eval.titv, Novelty == "all")$tiTvRatio;
	eval.titv.known = subset(eval.titv, Novelty == "known")$tiTvRatio;
	eval.titv.novel = subset(eval.titv, Novelty == "novel")$tiTvRatio;

	table4 = matrix(c(eval.counts.all, eval.counts.known, eval.counts.novel, eval.titv.all, eval.titv.known, eval.titv.novel, "3.0 - 3.2", "3.2 - 3.4", "2.7 - 3.0"), nrow=3);

	rownames(table4) = c("All", "Known", "Novel");
	colnames(table4) = c("Found", "Ti/Tv ratio", "Expected Ti/Tv ratio");

	

	par(mar=c(0,0,0,0))
	textplot(table4, rmar=1, col.rownames="dark blue", cex=1.25, valign="top")
	title(main="Variant Summary", family="sans", cex.main=1.25, line=-2)
# 	
# #plots
# #fix this reader
# 	eval.bysample = read.csv(paste(cmdargs$evalroot, ".SimpleMetricsBySample.csv", sep=""), header=TRUE, comment.char="#");
# 	eval.bysample.called = subset(eval.bysample, evaluation_name == "eval" & comparison_name == "dbsnp" & jexl_expression == "none" & filter_name == "called");
# 	eval.bysample.all = subset(eval.bysample.called, novelty_name == "all");
# 	eval.bysample.known = subset(eval.bysample.called, novelty_name == "known");
# 	eval.bysample.novel = subset(eval.bysample.called, novelty_name == "novel");

	eval.ac = eval$SimpleMetricsByAC.metrics
	eval.ac.all = subset(eval.ac, Novelty == "all");
	eval.ac.known = subset(eval.ac, Novelty == "known");
 	eval.ac.novel = subset(eval.ac, Novelty == "novel");
# 
# 	eval.func = read.csv(paste(cmdargs$evalroot, ".Functional_Class_Counts_by_Sample.csv", sep=""), header=TRUE, comment.char="#");
# 	eval.func.called = subset(eval.func, evaluation_name == "eval" & comparison_name == "dbsnp" & jexl_expression == "none" & filter_name == "called");
# 	eval.func.all = subset(eval.func.called, novelty_name == "all");
# 	eval.func.known = subset(eval.func.called, novelty_name == "known");
# 	eval.func.novel = subset(eval.func.called, novelty_name == "novel");


	#boxplot(eval.bysample.all$CountVariants, eval.bysample.known$CountVariants, eval.bysample.novel$CountVariants, names=c("All", "Known", "Novel"), ylab="Variants per sample", main="", cex=1.3, cex.lab=1.3, cex.axis=1.3);

# 	par(mar=c(5, 4, 4, 2) + 0.1)
# 	ind = order(eval.bysample.all$CountVariants);
# 	plot(c(1:length(eval.bysample.all$CountVariants)), eval.bysample.all$CountVariants[ind], col="black", cex=1.1, cex.lab=1.1, cex.axis=1.1, main="Variants per Sample", xlab="Sample", ylab="Number of variants", bty="n", ylim=c(0, max(eval.bysample.all$CountVariants)));
# 	points(c(1:length(eval.bysample.known$CountVariants)), eval.bysample.known$CountVariants[ind], col="blue", cex=1.3);
# 	points(c(1:length(eval.bysample.novel$CountVariants)), eval.bysample.novel$CountVariants[ind], col="red", cex=1.3);
# 	legend("right", max(eval.bysample.all$CountVariants)/2, c("All", "Known", "Novel"), col=c("black", "blue", "red"), pt.cex=1.3, pch=21);

	par(mar=c(5, 4, 4, 2) + 0.1)
	plot(eval.ac.all$AC, eval.ac.all$n, col="black", type="l", lwd=2, cex=1.1, cex.lab=1.1, cex.axis=1.1, xlab="Allele count", ylab="Number of variants", main="Variants by Allele Count", log="xy", bty="n");
	points(eval.ac.known$AC, eval.ac.known$n, col="blue", type="l", lwd=2);
	points(eval.ac.novel$AC, eval.ac.novel$n, col="red", type="l", lwd=2);
	legend("topright", c("All", "Known", "Novel"), col=c("black", "blue", "red"), lwd=2);

	#plot(eval.func.all$Synonymous[ind] / (eval.func.all$Missense + eval.func.all$Nonsense)[ind], ylim=c(0, 2), cex=1.3, cex.lab=1.3, cex.axis=1.3, bty="n", xlab="Sample", ylab="Ratio of synonymous to non-synonymous variants", col="black");
	#points(eval.func.known$Synonymous[ind] / (eval.func.known$Missense + eval.func.known$Nonsense)[ind], cex=1.3, col="blue");
	#points(eval.func.novel$Synonymous[ind] / (eval.func.novel$Missense + eval.func.novel$Nonsense)[ind], cex=1.3, col="red");
	#legend("topright", c("All", "Known", "Novel"), col=c("black", "blue", "red"), pt.cex=1.3, pch=21);


	
	dev.off()
	}

tearsheet()

# Plots 
plots<-function(){
# eval.bysample = read.csv(paste(cmdargs$evalroot, ".SimpleMetricsBySample.csv", sep=""), header=TRUE, comment.char="#");
# eval.bysample.called = subset(eval.bysample, evaluation_name == "eval" & comparison_name == "dbsnp" & jexl_expression == "none" & filter_name == "called");
# eval.bysample.all = subset(eval.bysample.called, novelty_name == "all");
# eval.bysample.known = subset(eval.bysample.called, novelty_name == "known");
# eval.bysample.novel = subset(eval.bysample.called, novelty_name == "novel");


  eval.ac = eval$SimpleMetricsByAC.metrics
	eval.ac.all = subset(eval.ac.called, Novelty == "all");
	eval.ac.known = subset(eval.ac.called, Novelty == "known");
 	eval.ac.novel = subset(eval.ac.called, Novelty == "novel");
# 
# eval.func = read.csv(paste(cmdargs$evalroot, ".Functional_Class_Counts_by_Sample.csv", sep=""), header=TRUE, comment.char="#");
# eval.func.called = subset(eval.func, evaluation_name == "eval" & comparison_name == "dbsnp" & jexl_expression == "none" & filter_name == "called");
# eval.func.all = subset(eval.func.called, novelty_name == "all");
# eval.func.known = subset(eval.func.called, novelty_name == "known");
# eval.func.novel = subset(eval.func.called, novelty_name == "novel");

  pdf(file= cmdargs$plotout, width=22, height=17, pagecentre=TRUE, pointsize=24) 
# 
# boxplot(eval.bysample.all$CountVariants, eval.bysample.known$CountVariants, eval.bysample.novel$CountVariants, names=c("All", "Known", "Novel"), ylab="Variants per sample", main="", cex=1.3, cex.lab=1.3, cex.axis=1.3);
# 
# ind = order(eval.bysample.all$CountVariants);
# plot(c(1:length(eval.bysample.all$CountVariants)), eval.bysample.all$CountVariants[ind], col="black", cex=1.3, cex.lab=1.3, cex.axis=1.3, xlab="Sample", ylab="Number of variants", bty="n", ylim=c(0, max(eval.bysample.all$CountVariants)));
# points(c(1:length(eval.bysample.known$CountVariants)), eval.bysample.known$CountVariants[ind], col="blue", cex=1.3);
# points(c(1:length(eval.bysample.novel$CountVariants)), eval.bysample.novel$CountVariants[ind], col="red", cex=1.3);
# legend(0, max(eval.bysample.all$CountVariants)/2, c("All", "Known", "Novel"), col=c("black", "blue", "red"), pt.cex=1.3, pch=21);

plot(eval.ac.all$AC, eval.ac.all$n, col="black", type="l", lwd=2, cex=1.3, cex.lab=1.3, cex.axis=1.3, xlab="Allele count", ylab="Number of variants", main="", log="xy", bty="n");
points(eval.ac.known$AC, eval.ac.known$n, col="blue", type="l", lwd=2);
points(eval.ac.novel$AC, eval.ac.novel$n, col="red", type="l", lwd=2);
legend("topright", c("All", "Known", "Novel"), col=c("black", "blue", "red"), lwd=2);
# 
# plot(eval.func.all$Synonymous[ind] / (eval.func.all$Missense + eval.func.all$Nonsense)[ind], ylim=c(0, 2), cex=1.3, cex.lab=1.3, cex.axis=1.3, bty="n", xlab="Sample", ylab="Ratio of synonymous to non-synonymous variants", col="black");
# points(eval.func.known$Synonymous[ind] / (eval.func.known$Missense + eval.func.known$Nonsense)[ind], cex=1.3, col="blue");
# points(eval.func.novel$Synonymous[ind] / (eval.func.novel$Missense + eval.func.novel$Nonsense)[ind], cex=1.3, col="red");
# legend("topright", c("All", "Known", "Novel"), col=c("black", "blue", "red"), pt.cex=1.3, pch=21);

dev.off();
}
