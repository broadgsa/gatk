args = commandArgs(TRUE);

evalRoot = args[1];
if (is.na(evalRoot)) { evalRoot = "results/v2/intermediate/CEU.combined.vcf.eval/eval"; }

plotRoot = args[2];
if (is.na(plotRoot)) { plotRoot = "results/v2/intermediate/CEU.combined.vcf.eval/eval.plot"; }

population = args[3];
if (is.na(population)) { population = "unknown population"; }

fileAlleleCountStats = paste(evalRoot, ".AlleleCountStats.csv", sep="");
fileCompOverlap = paste(evalRoot, ".Comp_Overlap.csv", sep="");
fileCountVariants = paste(evalRoot, ".Count_Variants.csv", sep="");
fileGenotypeConcordance = paste(evalRoot, ".Genotype_Concordance.csv", sep="");
fileMetricsByAc = paste(evalRoot, ".MetricsByAc.csv", sep="");
fileQuality_Metrics_by_allele_count = paste(evalRoot, ".Quality_Metrics_by_allele_count.csv", sep="");
fileQualityScoreHistogram = paste(evalRoot, ".QualityScoreHistogram.csv", sep="");
fileSampleStatistics = paste(evalRoot, ".Sample_Statistics.csv", sep="");
fileSampleSummaryStatistics = paste(evalRoot, ".Sample_Summary_Statistics.csv", sep="");
fileTi_slash_Tv_Variant_Evaluator = paste(evalRoot, ".Ti_slash_Tv_Variant_Evaluator.csv", sep="");
fileTiTvStats = paste(evalRoot, ".TiTvStats.csv", sep="");
fileVariant_Quality_Score = paste(evalRoot, ".Variant_Quality_Score.csv", sep="");
filePerSample = paste(evalRoot, ".persample.table", sep="");

# Allele count spectrum
if (!is.na(plotRoot)) {
    pdf(paste(plotRoot, ".acs.pdf", sep=""));
} else {
    x11();
}

dataMetricsByAc = read.csv(fileMetricsByAc, header=TRUE, comment.char="#");

dataMetricsByAc.none.all = dataMetricsByAc[which(dataMetricsByAc$jexl_expression == "none" & dataMetricsByAc$filter_name == "called" & dataMetricsByAc$novelty_name == "all"),];

dataMetricsByAc.int.all = dataMetricsByAc[which(dataMetricsByAc$jexl_expression == "Intersection" & dataMetricsByAc$filter_name == "called" & dataMetricsByAc$novelty_name == "all"),];
dataMetricsByAc.int.known = dataMetricsByAc[which(dataMetricsByAc$jexl_expression == "Intersection" & dataMetricsByAc$filter_name == "called" & dataMetricsByAc$novelty_name == "known"),];
dataMetricsByAc.int.novel = dataMetricsByAc[which(dataMetricsByAc$jexl_expression == "Intersection" & dataMetricsByAc$filter_name == "called" & dataMetricsByAc$novelty_name == "novel"),];

dataMetricsByAc.UG.all = dataMetricsByAc[which(dataMetricsByAc$jexl_expression == "UG" & dataMetricsByAc$filter_name == "called" & dataMetricsByAc$novelty_name == "all"),] + dataMetricsByAc.int.all;
dataMetricsByAc.UG.known = dataMetricsByAc[which(dataMetricsByAc$jexl_expression == "UG" & dataMetricsByAc$filter_name == "called" & dataMetricsByAc$novelty_name == "known"),] + dataMetricsByAc.int.known;
dataMetricsByAc.UG.novel = dataMetricsByAc[which(dataMetricsByAc$jexl_expression == "UG" & dataMetricsByAc$filter_name == "called" & dataMetricsByAc$novelty_name == "novel"),] + dataMetricsByAc.int.novel;

dataMetricsByAc.QCALL.all = dataMetricsByAc[which(dataMetricsByAc$jexl_expression == "QCALL" & dataMetricsByAc$filter_name == "called" & dataMetricsByAc$novelty_name == "all"),] + dataMetricsByAc.int.all;
dataMetricsByAc.QCALL.known = dataMetricsByAc[which(dataMetricsByAc$jexl_expression == "QCALL" & dataMetricsByAc$filter_name == "called" & dataMetricsByAc$novelty_name == "known"),] + dataMetricsByAc.int.known;
dataMetricsByAc.QCALL.novel = dataMetricsByAc[which(dataMetricsByAc$jexl_expression == "QCALL" & dataMetricsByAc$filter_name == "called" & dataMetricsByAc$novelty_name == "novel"),] + dataMetricsByAc.int.novel;


plot(0, 0, type="n", xlab="Allele count", ylab="Number of variants", xlim=c(0, max(dataMetricsByAc.none.all$AC)), ylim=c(min(dataMetricsByAc.none.all$n), max(dataMetricsByAc.none.all$n)), main=paste("Allele count spectrum (", population, ")", sep=""), cex=1.2, cex.axis=1.2, cex.lab=1.2, bty="n");

points(dataMetricsByAc.int.all$AC, dataMetricsByAc.int.all$n, type="l", lwd=3, col="black");
points(dataMetricsByAc.int.known$AC, dataMetricsByAc.int.known$n, type="l", lwd=3, col="blue");
points(dataMetricsByAc.int.novel$AC, dataMetricsByAc.int.novel$n, type="l", lwd=3, col="red");

points(dataMetricsByAc.UG.all$AC, dataMetricsByAc.UG.all$n, type="l", lwd=1, lty=2, col="black");
points(dataMetricsByAc.UG.known$AC, dataMetricsByAc.UG.known$n, type="l", lwd=1, lty=2, col="blue");
points(dataMetricsByAc.UG.novel$AC, dataMetricsByAc.UG.novel$n, type="l", lwd=1, lty=2, col="red");

points(dataMetricsByAc.QCALL.all$AC, dataMetricsByAc.QCALL.all$n, type="l", lwd=1, lty=3, col="black");
points(dataMetricsByAc.QCALL.known$AC, dataMetricsByAc.QCALL.known$n, type="l", lwd=1, lty=3, col="blue");
points(dataMetricsByAc.QCALL.novel$AC, dataMetricsByAc.QCALL.novel$n, type="l", lwd=1, lty=3, col="red");

legend("topright", c("Intersection: all", "Intersection: known (present in dbSNP 129)", "Intersection: novel (absent in dbSNP 129)", "UnifiedGenotyper: all", "UnifiedGenotyper: known (present in dbSNP 129)", "UnifiedGenotyper: novel (absent in dbSNP 129)", "QCALL: all", "QCALL: known (present in dbSNP 129)", "QCALL: novel (absent in dbSNP 129)"), lwd=c(3, 3, 3, 1, 1, 1, 1, 1, 1), lty=c(1, 1, 1, 2, 2, 2, 3, 3, 3), col=c("black", "blue", "red"), cex=1.2);

if (!is.na(plotRoot)) { dev.off(); }

# Ti/Tv per allele count
if (!is.na(plotRoot)) {
    pdf(paste(plotRoot, ".titv.pdf", sep=""));
} else {
    x11();
}

plot(0, 0, type="n", xlab="Allele count", ylab="Transition/transversion ratio", xlim=c(0, max(dataMetricsByAc.none.all$AC)), ylim=c(0, max(dataMetricsByAc$Ti.Tv)), main=paste("Ti/Tv per allele count (", population, ")", sep=""), cex=1.2, cex.lab=1.2, cex.axis=1.2, bty="n");

points(dataMetricsByAc.int.all$AC, dataMetricsByAc.int.all$Ti.Tv, type="l", lwd=3, lty=1, col="black");
points(dataMetricsByAc.int.known$AC, dataMetricsByAc.int.known$Ti.Tv, type="l", lwd=3, lty=1, col="blue");
points(dataMetricsByAc.int.novel$AC, dataMetricsByAc.int.novel$Ti.Tv, type="l", lwd=3, lty=1, col="red");

points(dataMetricsByAc.UG.all$AC, dataMetricsByAc.int.all$Ti.Tv, type="l", lwd=1, lty=2, col="black");
points(dataMetricsByAc.UG.known$AC, dataMetricsByAc.int.known$Ti.Tv, type="l", lwd=1, lty=2, col="blue");
points(dataMetricsByAc.UG.novel$AC, dataMetricsByAc.int.novel$Ti.Tv, type="l", lwd=1, lty=2, col="red");

points(dataMetricsByAc.QCALL.all$AC, dataMetricsByAc.int.all$Ti.Tv, type="l", lwd=1, lty=3, col="black");
points(dataMetricsByAc.QCALL.known$AC, dataMetricsByAc.int.known$Ti.Tv, type="l", lwd=1, lty=3, col="blue");
points(dataMetricsByAc.QCALL.novel$AC, dataMetricsByAc.int.novel$Ti.Tv, type="l", lwd=1, lty=3, col="red");

legend("topright", c("Intersection: all", "Intersection: known (present in dbSNP 129)", "Intersection: novel (absent in dbSNP 129)", "UnifiedGenotyper: all", "UnifiedGenotyper: known (present in dbSNP 129)", "UnifiedGenotyper: novel (absent in dbSNP 129)", "QCALL: all", "QCALL: known (present in dbSNP 129)", "QCALL: novel (absent in dbSNP 129)"), lwd=c(3, 3, 3, 1, 1, 1, 1, 1, 1), lty=c(1, 1, 1, 2, 2, 2, 3, 3, 3), col=c("black", "blue", "red"), cex=1.2);

if (!is.na(plotRoot)) { dev.off(); }

# Barplot (unnormalized)
if (!is.na(plotRoot)) {
    pdf(paste(plotRoot, ".venn_per_ac.pdf", sep=""));
} else {
    x11();
}

norm = dataMetricsByAc.UG.all$n + dataMetricsByAc.int.all$n + dataMetricsByAc.QCALL.all$n;
mat = rbind(dataMetricsByAc.QCALL.all$n, dataMetricsByAc.int.all$n, dataMetricsByAc.UG.all$n);
matnorm = rbind(dataMetricsByAc.QCALL.all$n/norm, dataMetricsByAc.int.all$n/norm, dataMetricsByAc.UG.all$n/norm);

qcallcolor = "#FF5555";
intcolor = "#5582C6";
ugcolor = "#55BBFF";

barplot(mat, col=c(qcallcolor, intcolor, ugcolor), xlab="Allele count", ylab="counts", main=paste("Callset concordance per allele count (", population, ")", sep=""), names.arg=dataMetricsByAc.UG.all$AC);
legend("topright", c("UG-only", "Intersection", "QCALL-only"), fill=c(ugcolor, intcolor, qcallcolor));

if (!is.na(plotRoot)) { dev.off(); }

# Barplot (normalized)
if (!is.na(plotRoot)) {
    pdf(paste(plotRoot, ".venn_per_ac.norm.pdf", sep=""));
} else {
    x11();
}

barplot(matnorm, col=c(qcallcolor, intcolor, ugcolor), xlab="Allele count", ylab="Fraction", main=paste("Callset concordance per allele count (", population, ")", sep=""), names.arg=dataMetricsByAc.UG.all$AC, ylim=c(0, 1.3));
legend("topright", c("UG-only", "Intersection", "QCALL-only"), fill=c(ugcolor, intcolor, qcallcolor));

if (!is.na(plotRoot)) { dev.off(); }

# Per-sample
if (!is.na(plotRoot)) {
    pdf(paste(plotRoot, ".counts_per_sample.pdf", sep=""), width=10);
} else {
    x11(width=10);
}

dataPerSample = read.table(filePerSample, header=TRUE);

par(mar=c(7, 5, 2, 1));

plot(0, 0, type="n", xaxt="n", xlim=c(1, length(dataPerSample$sample)), ylim=c(0, max(dataPerSample$Intersection + dataPerSample$UG, dataPerSample$Intersection, dataPerSample$Intersection + dataPerSample$QCALL)), xlab="", ylab="Counts", main=paste("SNPs per sample (", population, ")", sep=""), cex=1.2, cex.lab=1.2, cex.axis=1.2);

points(dataPerSample$sample, dataPerSample$Intersection + dataPerSample$UG, col=ugcolor);
points(dataPerSample$sample, dataPerSample$Intersection, col=intcolor);
points(dataPerSample$sample, dataPerSample$Intersection + dataPerSample$QCALL, col=qcallcolor);

legend("bottomleft", c("UG", "Intersection", "QCALL"), fill=c(ugcolor, intcolor, qcallcolor));

axis(1, at=c(1:length(dataPerSample$sample)), labels=dataPerSample$sample, las=2, cex.axis=0.7);

if (!is.na(plotRoot)) { dev.off(); }
