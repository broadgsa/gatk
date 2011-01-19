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
        bamlist = list(value=NA, doc="list of BAM files"),
        evalroot = list(value=NA, doc="VariantEval root"),
        tearout = list(value=NA, doc="Output path for tearsheet PDF"),
        plotout = list(value=NA, doc="Output path for PDF")
    ),
    doc="Creates a tearsheet"
);

#if (0) {
bamlist = read.table(cmdargs$bamlist);

fclanes = c();
for (bam in bamlist$V1) {
    bamheader = system(paste("samtools view -H", bam), intern=TRUE);

    if (length(bamheader) > 0) {
        rgs = bamheader[grep("^@RG", bamheader)];

        for (rg in rgs) {
            id = grep("ID:", unlist(strsplit(rg, "\t")), value=TRUE);
            id = sub("ID:", "", id);

            fclanes = c(fclanes, id);
        }
    } else {
        Print(sprintf("Could not load '%s'\n", bam));
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
squid_fclanes = gsub("A.XX", "", squid_fclanes);

dproj = d[which(squid_fclanes %in% fclanes),];
d2proj = d2[which(d2$"Project" %in% unique(dproj$"Project") & d2$"Sample" %in% dproj$"External ID"),];
#}

tearsheetdrop <- "tearsheetdrop.jpg" #put the path to the tearsheet backdrop here

tearsheet<-function(){
	
	pdf(file= cmdargs$tearout, width=22, height=15, pagecentre=TRUE, pointsize=24) 
	
	#define layout
	layout(matrix(c(1,1,2,4,3,5), ncol=2, nrow=3, byrow=TRUE), heights=c(1, 2.5,2.5), respect=FALSE)
    
    #prep for title bar
    full<-strsplit(cmdargs$evalroot, "/")
    name<-strsplit(full[[1]][length(full[[1]])], ".",fixed=TRUE)[[1]][1]
    title=paste(name, ": TEAR SHEET", sep="")
    drop<-read.jpeg(system.file(tearsheetdrop, package="gsalib"))
    
    #plot title bar
    par(mar=c(0,0,0,0))
    plot(drop)
    text(110, 45, title, family="serif", adj=c(0,0), cex=3, col=gray(.25))
    

	# Project summary
	projects = paste(unique(dproj$"Project"), collapse=", ");

	used_samples = length(bamlist$V1);

	unused_samples = 0;

	sequencing_protocol = "Hybrid selection"; #can this be extracted?

	bait_design = paste(unique(dproj$"Bait Set"), collapse=", ");

	callable_target = paste(unique(dproj$"Target Territory"), collapse=", ");

	table1<-rbind(paste(used_samples," used samples/", unused_samples + used_samples," total samples", sep=""), sequencing_protocol, bait_design, callable_target)
	print(nrow(table1))
	rownames(table1)<-c("Samples","Sequencing Protocol", "Bait Design","Callable Target")
	par(mar=c(4,4,4,4))
	textplot(table1, col.rownames="darkblue", show.colnames=FALSE, cex=1.25)
    title(main=sprintf("Project Summary (%s)\n", projects), family="sans", cex.main=1.75)

	# Bases summary 

	reads_per_lane_mean = signif(mean(dproj$"PF Reads (HS)"), 3);
	reads_per_lane_sd = signif(sd(dproj$"PF Reads (HS)"), 3);
	lanes<-sprintf("%s +/- %s\n", reads_per_lane_mean, reads_per_lane_sd)

	used_bases_per_lane_mean = signif(mean(dproj$"PF HQ Aligned Q20 Bases"), 3);
	used_bases_per_lane_sd = signif(sd(dproj$"PF HQ Aligned Q20 Bases"), 3);
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


	reads_per_sample_mean = signif(mean(d2proj$"PF Reads"), 3);
	reads_per_sample_sd = signif(sd(d2proj$"PF Reads"), 3);
	samps<-sprintf("%s +/- %s\n", reads_per_sample_mean, reads_per_sample_sd);

	used_bases_per_sample_mean = signif(mean(d2proj$"PF HQ Aligned Q20 Bases"), 3);
	used_bases_per_sample_sd = signif(sd(d2proj$"PF HQ Aligned Q20 Bases"), 3);
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
		print(nrow(table2))

	rownames(table2)<-c("Reads", "Used bases", "Average target coverage", "% loci covered to 10x", "% loci covered to 20x","% loci covered to 30x")
	par(mar=c(4,4,4,4))
	textplot(table2, rmar=1, col.rownames="dark blue", cex=1.25)
	title(main="Bases Summary", family="sans", cex.main=1.75)


# Sequencing summary

	instrument <- c();
	if(length(grep("AAXX", dproj$Flowcell))>0){
		instrument <- c(instrument, "Illumina GA2")
		}
	if(length(grep("ABXX", dproj$Flowcell))>0){
		instrument <- c(instrument, "Illumina HiSeq")
		}

	used_lanes = nrow(dproj);
	unused_lanes_by_sequencing = 0; #can we get this?
	unused_lanes_by_analysis = 0;


	lanes_per_sample_mean = mean(table(dproj$"External ID"));
	lanes_per_sample_sd = sd(table(dproj$"External ID"));
	lanes_per_sample_median = median(table(dproj$"External ID"));
	lanes_paired = nrow(subset(dproj, dproj$"Lane Type" == "Paired"));
	lanes_widowed = nrow(subset(dproj, dproj$"Lane Type" == "Widowed"));
	lanes_single = nrow(subset(dproj, dproj$"Lane Type" == "Single"));

	read_length_mean = mean(dproj$"Mean Read Length (P)");
	read_length_sd = sd(dproj$"Mean Read Length (P)");
	read_length_median = median(dproj$"Mean Read Length (P)");

	date = dproj$"Run Date";
	date = sub("JAN", "01", date);
	date = sub("FEB", "02", date);
	date = sub("MAR", "03", date);
	date = sub("APR", "04", date);
	date = sub("MAY", "05", date);
	date = sub("JUN", "06", date);
	date = sub("JUL", "07", date);
	date = sub("AUG", "08", date);
	date = sub("SEP", "09", date);
	date = sub("OCT", "10", date);
	date = sub("NOV", "11", date);
	date = sub("DEC", "12", date);
	date = date[order(as.Date(date, format="%d-%m-%Y"))];

	start_date = date[1];
	end_date = date[length(date)];

	
	table3<-rbind(paste(instrument), used_lanes, sprintf("%s rejected by sequencing, %s by analysis\n", unused_lanes_by_sequencing, unused_lanes_by_analysis), sprintf("%0.1f +/- %0.1f lanes (median=%0.1f)\n", lanes_per_sample_mean, lanes_per_sample_sd, lanes_per_sample_median), sprintf("%s paired, %s widowed, %s single\n", lanes_paired, lanes_widowed, lanes_single), sprintf("%0.1f +/- %0.1f bases (median=%0.1f)\n", read_length_mean, read_length_sd, read_length_median), sprintf("\tSequencing dates: %s to %s\n", start_date, end_date))
		print(nrow(table3))

  	rownames(table3)<-c("Sequencer", "Used lanes", "Unused lanes","Used lanes/sample", "Lane pariteies", "Read legnths", "Sequencing dates")
	par(mar=c(4,4,4,4))
	textplot(table3, rmar=1, col.rownames="dark blue", show.colnames=FALSE, cex=1.25)
	title(main="Sequencing Summary", family="sans", cex.main=1.75)

# Variant summary
	
	eval.counts = read.csv(paste(cmdargs$evalroot, ".Count_Variants.csv", sep=""), header=TRUE, comment.char="#");
	eval.counts.called = subset(eval.counts, evaluation_name == "eval" & comparison_name == "dbsnp" & 	jexl_expression == "none" & filter_name == "called");
	eval.counts.called.all = subset(eval.counts.called, novelty_name == "all")$nVariantLoci;
	eval.counts.called.known = subset(eval.counts.called, novelty_name == "known")$nVariantLoci;
	eval.counts.called.novel = subset(eval.counts.called, novelty_name == "novel")$nVariantLoci;

	eval.titv = read.csv(paste(cmdargs$evalroot, ".Ti_slash_Tv_Variant_Evaluator.csv", sep=""), header=TRUE, comment.char="#");
	eval.titv.called = subset(eval.titv, evaluation_name == "eval" & comparison_name == "dbsnp" & 	jexl_expression == "none" & filter_name == "called");
	eval.titv.called.all = subset(eval.titv.called, novelty_name == "all")$ti.tv_ratio;
	eval.titv.called.known = subset(eval.titv.called, novelty_name == "known")$ti.tv_ratio;
	eval.titv.called.novel = subset(eval.titv.called, novelty_name == "novel")$ti.tv_ratio;

	table4 = matrix(c(eval.counts.called.all, eval.counts.called.known, eval.counts.called.novel, eval.titv.called.all, eval.titv.called.known, eval.titv.called.novel, "3.0 - 3.2", "3.2 - 3.4", "2.7 - 3.0"), nrow=3);
		print(nrow(table4))

	rownames(table4) = c("All", "Known", "Novel");
	colnames(table4) = c("Found", "Ti/Tv ratio", "Expected Ti/Tv ratio");

	

	textplot(table4, rmar=1, col.rownames="dark blue", cex=1.25)
	title(main="Variant Summary", family="sans", cex.main=1.75)
	
	dev.off()
	}

tearsheet()

# Plots 
plots<-function(){
eval.bysample = read.csv(paste(cmdargs$evalroot, ".SimpleMetricsBySample.csv", sep=""), header=TRUE, comment.char="#");
eval.bysample.called = subset(eval.bysample, evaluation_name == "eval" & comparison_name == "dbsnp" & jexl_expression == "none" & filter_name == "called");
eval.bysample.called.all = subset(eval.bysample.called, novelty_name == "all");
eval.bysample.called.known = subset(eval.bysample.called, novelty_name == "known");
eval.bysample.called.novel = subset(eval.bysample.called, novelty_name == "novel");

eval.ac = read.csv(paste(cmdargs$evalroot, ".MetricsByAc.csv", sep=""), header=TRUE, comment.char="#");
eval.ac.called = subset(eval.ac, evaluation_name == "eval" & comparison_name == "dbsnp" & jexl_expression == "none" & filter_name == "called");
eval.ac.called.all = subset(eval.ac.called, novelty_name == "all");
eval.ac.called.known = subset(eval.ac.called, novelty_name == "known");
eval.ac.called.novel = subset(eval.ac.called, novelty_name == "novel");

eval.func = read.csv(paste(cmdargs$evalroot, ".Functional_Class_Counts_by_Sample.csv", sep=""), header=TRUE, comment.char="#");
eval.func.called = subset(eval.func, evaluation_name == "eval" & comparison_name == "dbsnp" & jexl_expression == "none" & filter_name == "called");
eval.func.called.all = subset(eval.func.called, novelty_name == "all");
eval.func.called.known = subset(eval.func.called, novelty_name == "known");
eval.func.called.novel = subset(eval.func.called, novelty_name == "novel");

pdf(cmdargs$plotout);

boxplot(eval.bysample.called.all$CountVariants, eval.bysample.called.known$CountVariants, eval.bysample.called.novel$CountVariants, names=c("All", "Known", "Novel"), ylab="Variants per sample", main="", cex=1.3, cex.lab=1.3, cex.axis=1.3);

ind = order(eval.bysample.called.all$CountVariants);
plot(c(1:length(eval.bysample.called.all$CountVariants)), eval.bysample.called.all$CountVariants[ind], col="black", cex=1.3, cex.lab=1.3, cex.axis=1.3, xlab="Sample", ylab="Number of variants", bty="n", ylim=c(0, max(eval.bysample.called.all$CountVariants)));
points(c(1:length(eval.bysample.called.known$CountVariants)), eval.bysample.called.known$CountVariants[ind], col="blue", cex=1.3);
points(c(1:length(eval.bysample.called.novel$CountVariants)), eval.bysample.called.novel$CountVariants[ind], col="red", cex=1.3);
legend(0, max(eval.bysample.called.all$CountVariants)/2, c("All", "Known", "Novel"), col=c("black", "blue", "red"), pt.cex=1.3, pch=21);

plot(eval.ac.called.all$AC, eval.ac.called.all$n, col="black", type="l", lwd=2, cex=1.3, cex.lab=1.3, cex.axis=1.3, xlab="Allele count", ylab="Number of variants", main="", log="xy", bty="n");
points(eval.ac.called.known$AC, eval.ac.called.known$n, col="blue", type="l", lwd=2);
points(eval.ac.called.novel$AC, eval.ac.called.novel$n, col="red", type="l", lwd=2);
legend("topright", c("All", "Known", "Novel"), col=c("black", "blue", "red"), lwd=2);

plot(eval.func.called.all$Synonymous[ind] / (eval.func.called.all$Missense + eval.func.called.all$Nonsense)[ind], ylim=c(0, 2), cex=1.3, cex.lab=1.3, cex.axis=1.3, bty="n", xlab="Sample", ylab="Ratio of synonymous to non-synonymous variants", col="black");
points(eval.func.called.known$Synonymous[ind] / (eval.func.called.known$Missense + eval.func.called.known$Nonsense)[ind], cex=1.3, col="blue");
points(eval.func.called.novel$Synonymous[ind] / (eval.func.called.novel$Missense + eval.func.called.novel$Nonsense)[ind], cex=1.3, col="red");
legend("topright", c("All", "Known", "Novel"), col=c("black", "blue", "red"), pt.cex=1.3, pch=21);

dev.off();
}

plots()
