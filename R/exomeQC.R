library("gsalib", lib.loc="/Users/depristo/Desktop/broadLocal/GATK/trunk/R/")
require("ggplot2")
require("gplots")

# TODOs:
#  Assumes you have indels in your call set.  If not you will get errors
#  Create pre/post calling sections
#  Allow conditional use of the preQCFile (where it's not available)

args = commandArgs(TRUE)
onCMDLine = ! is.na(args[1])
LOAD_DATA = T

# creates an array of c(sampleName1, ..., sampleNameN)
parseHighlightSamples <- function(s) {
  return(unlist(strsplit(s, ",", fixed=T)))
}

preQCFile = NA
if ( onCMDLine ) {
  ProjectName = args[1]
  VariantEvalRoot = args[2]
  outputPDF = args[3]
  if ( ! is.na(args[4]) ) 
    preQCFile = args[4]
  if ( ! is.na(args[5]) ) 
    highlightSamples = parseHighlightSamples(args[5])
  else
    highlightSamples = c()
} else {
  ProjectName = "InDevelopmentInR"
  preQCFile <- NA # "~/Desktop/broadLocal/GATK/trunk/qcTestData/GoT2D_exomes_batch_005_per_sample_metrics.tsv"
  #VariantEvalRoot <- "qcTestData//ESPGO_Gabriel_NHLBI_eomi_june_2011_batch1"
  VariantEvalRoot <- "qcTestData/MC_Engle_11_Samples_06092011"
  outputPDF = "bar.pdf"
  highlightSamples = c() # parseHighlightSamples("29029,47243")
}

print("Report")
print(paste("Project          :", ProjectName))
print(paste("VariantEvalRoot  :", VariantEvalRoot))
print(paste("outputPDF        :", outputPDF))
print(paste("preQCFile        :", preQCFile))
print(paste("highlightSamples :", highlightSamples))

expandVEReport <- function(d) {
  d$TiTvVariantEvaluator$tiTvRatio = round(d$TiTvVariantEvaluator$tiTvRatio,2) 
  d$CountVariants$deletionInsertionRatio = round(d$CountVariants$deletionInsertionRatio,2) 
  d$CountVariants$nIndels = d$CountVariants$nInsertions + d$CountVariants$nDeletions
  return(d)
}

# -------------------------------------------------------
# Utilities for displaying multiple plots per page
# -------------------------------------------------------

# Viewport (layout 2 graphs top to bottom)
distributePerSampleGraph <- function(distgraph, perSampleGraph, heights = c(2,1)) {
  Layout <- grid.layout(nrow = 2, ncol = 1, heights=heights)
  grid.newpage()
  pushViewport(viewport(layout = Layout))
  subplot <- function(x) viewport(layout.pos.row = x, layout.pos.col = 1)
  print(perSampleGraph, vp = subplot(1))
  print(distgraph, vp = subplot(2))
}

createMetricsBySites <- function(VariantEvalRoot, PreQCMetrics) {
  # Metrics by sites:
  #  bySite -> counts of SNPs and Indels by novelty, with expectations
  #  byAC -> snps and indels (known / novel)
  r = list( bySite = expandVEReport(gsa.read.gatkreport(paste(VariantEvalRoot, ".summary.eval", sep=""))),
               byAC = gsa.read.gatkreport(paste(VariantEvalRoot, ".byAC.eval", sep="")))
  r$byAC$CountVariants$nIndels = r$byAC$CountVariants$nInsertions + r$byAC$CountVariants$nDeletions
  r$byAC$TiTvVariantEvaluator$nSNPs = r$byAC$TiTvVariantEvaluator$nTi + r$byAC$TiTvVariantEvaluator$nTv
  r$byAC$CountVariants$AC = r$byAC$CountVariants$AlleleCount
  r$byAC$TiTvVariantEvaluator$AC = r$byAC$TiTvVariantEvaluator$AlleleCount
  return(r)
}

summaryTable <- function(metricsBySites, metricsBySample) {
  # SNP summary statistics
  merged = merge(metricsBySites$bySite$CountVariants, metricsBySites$bySite$TiTvVariantEvaluator)
  sub <- subset(merged, FunctionalClass=="all")
  raw = melt(sub, id.vars=c("Novelty"), measure.vars=c("nProcessedLoci", "nSNPs", "tiTvRatio", "nIndels", "deletionInsertionRatio"))
  table = cast(raw, Novelty ~ ...)
  # doesn't work with textplot
  colnames(table) <- c("Novelty", "Target size (bp)", "No. SNPs", "Ti/Tv", "No. Indels", "deletion/insertion ratio")
  return(table)
}

sampleSummaryTable <- function(metricsBySample) {
  # SNP summary statistics
  raw <- melt(metricsBySamples, id.vars=c("Novelty", "Sample"), measure.vars=c("nProcessedLoci", "nSNPs", "tiTvRatio", "nIndels", "deletionInsertionRatio"))
  table = cast(raw, Novelty ~ variable, mean)
  table$nSNPs <- round(table$nSNPs, 0)
  table$nIndels <- round(table$nIndels, 0)
  table$tiTvRatio <- round(table$tiTvRatio, 2)
  table$deletionInsertionRatio <- round(table$deletionInsertionRatio, 2)
  colnames(table) <- c("Novelty", "Target size (bp)", "No. SNPs", "Ti/Tv", "No. Indels", "deletion/insertion ratio")
  return(table)
}

overallSummaryTable <- function(metricsBySites, metricsBySamples) {
  sitesSummary <- as.data.frame(summaryTable(metricsBySites, metricsBySamples))
  sitesSummary$Metric.Type <- "Sites"
  sampleSummary <- as.data.frame(sampleSummaryTable(metricsBySamples))
  sampleSummary$Metric.Type <- "Per-sample avg."
  # that last item puts the metric.type second in the list
  return(rbind(sitesSummary, sampleSummary)[, c(1,7,2,3,4,5,6)])
}

summaryPlots <- function(metricsBySites) {
  name = "SNP and Indel count by novelty and allele frequency" 
  molten = melt(subset(metricsBySites$byAC$CountVariants, Novelty != "all" & AC > 0), id.vars=c("Novelty", "AC"), measure.vars=c(c("nSNPs", "nIndels")))
  p <- ggplot(data=molten, aes(x=AC, y=value+1, color=Novelty, fill=Novelty), group=variable)
  p <- p + opts(title = name)
  p <- p + scale_y_log10("Number of variants")
  p <- p + geom_point(alpha=0.5, size=3)
  p <- p + geom_line(size=1)
  p <- p + facet_grid(variable ~ ., scales="free")
  p <- p + scale_x_continuous("Allele count (AC)")
  p2 <- p + scale_x_log10("Allele count (AC)")
  p2 <- p2 + opts(title = "")
  distributePerSampleGraph(p2, p, c(1,1))

  # Counts vs. Allele frequency 
  name = "Variant counts by allele count"
  for ( measure in c("nSNPs", "nIndels")) {
    molten = melt(subset(metricsBySites$byAC$CountVariants, AC > 0), id.vars=c("Novelty", "AC"), measure.vars=c(measure))
    p <- ggplot(data=molten, aes(x=AC, y=value+1, color=Novelty), group=variable)
    p <- p + opts(title = paste(name, ":", measure))
    p <- p + scale_y_log10("Number of variants")
    p <- p + scale_x_log10("Allele count (AC)")
    p <- p + geom_point(alpha=0.5, size=4)
    p <- p + geom_smooth(aes(weight=value), size=1, method="lm", formula = y ~ x)
    p <- p + facet_grid(Novelty ~ ., scales="free")
    print(p)
  }
  
  name = "Transition / transversion ratio by allele count"
  # nSNPs > 0 => requires that we have some data here, otherwise Ti/Tv is zero from VE  
  minSNPsToInclude = 0
  byACNoAll = subset(metricsBySites$byAC$TiTvVariantEvaluator, Novelty != "all" & AC > 0 & nSNPs > minSNPsToInclude)
  p <- ggplot(data=byACNoAll, aes(x=AC, y=tiTvRatio, color=Novelty))
  p <- p + scale_y_continuous("Transition / transversion ratio", limits=c(0,4))
  p <- p + opts(title = name)
  p <- p + geom_smooth(size=2)
  p <- p + geom_point(aes(size=log10(nSNPs), weight=nSNPs), alpha=0.5)
  p <- p + scale_x_continuous("Allele count (AC)")
  p2 <- p + scale_x_log10("Allele count (AC)")
  p2 <- p2 + opts(title = "")
  distributePerSampleGraph(p2, p, c(1,1))
  
  # SNPs to indels ratio by allele frequency
  name = "SNPs to indels ratio by allele frequency" 
  metricsBySites$byAC$CountVariants$SNP.Indel.Ratio = metricsBySites$byAC$CountVariants$nSNPs / metricsBySites$byAC$CountVariants$nIndels
  metricsBySites$byAC$CountVariants$SNP.Indel.Ratio[metricsBySites$byAC$CountVariants$nIndels == 0] = NaN
  p <- ggplot(data=subset(metricsBySites$byAC$CountVariants, Novelty == "all" & nSNPs > 0), aes(x=AC, y=SNP.Indel.Ratio))
  p <- p + opts(title = name)
  p <- p + scale_y_continuous("SNP to indel ratio")
  #p <- p + scale_y_log10()
  p <- p + geom_point(alpha=0.5, aes(size=log10(nIndels)))
  p <- p + geom_smooth(size=2, aes(weight=nIndels))
  print(p)
  
  name = "SNP counts by functional class" 
  molten = melt(subset(metricsBySites$bySite$CountVariants, Novelty != "all" & FunctionalClass != "all"), id.vars=c("Novelty", "FunctionalClass"), measure.vars=c(c("nSNPs")))
  p <- ggplot(data=molten, aes(x=FunctionalClass, y=value, fill=Novelty), group=FunctionalClass)
  p <- p + opts(title = name)
  p <- p + scale_y_log10("No. of SNPs")
  p <- p + geom_bar(position="dodge")
  print(p)
}

addSection <- function(name) {
    par("mar", c(5, 4, 4, 2))
    frame()
    title(name, cex=2)
}
 
# -------------------------------------------------------
# read functions
# -------------------------------------------------------

createMetricsBySamples <- function(VariantEvalRoot) {
  bySampleEval <- expandVEReport(gsa.read.gatkreport(paste(VariantEvalRoot, ".bySample.eval", sep="")))
  r = merge(bySampleEval$TiTvVariantEvaluator, bySampleEval$CountVariants)
  r = merge(r, bySampleEval$CompOverlap)
  if ( ! is.na(preQCFile) ) {
    preQCMetrics <- read.table(preQCFile, header=T)
    r = merge(r, preQCMetrics)
  }
  # order the samples by nSNPs -- it's the natural ordering.
  x = subset(r, Novelty=="all")
  r$Sample <- factor(x$Sample, levels=x$Sample[order(x$nSNPs)])

  # add highlight info
  r$highlight = r$Sample %in% highlightSamples

  #r = merge(merge(preQCMetrics, byACEval$TiTvVariantEvaluator), byACEval$CountVariants)
  return(subset(r, Sample != "all"))
}

# -------------------------------------------------------
# Per sample plots
# -------------------------------------------------------

perSamplePlots <- function(metricsBySamples) {
  metricsBySamples$highlightTextSizes = c(1,2)[metricsBySamples$highlight+1]
  sampleTextLabel <- geom_text(aes(label=Sample, size=highlightTextSizes)) 
  sampleTextLabelScale <- scale_size("Highlighted samples", to=c(3,5), breaks=c(1,2), labels=c("regular", "highlighted"))
  xAxis <- scale_x_discrete("Sample (ordered by nSNPs)", formatter=function(x) "")

  measures = c("nSNPs", "tiTvRatio", "nSingletons", "nIndels", "deletionInsertionRatio")
  name = "by sample"
  for ( measure in measures ) {
    molten = melt(metricsBySamples, id.vars=c("Novelty", "Sample", "highlightTextSizes"), measure.vars=c(measure))

    # distribution
    p1 <- ggplot(data=molten, aes(x=value, group=Novelty, fill=Novelty))
    #p1 <- p1 + opts(title = paste(measure, name))
    p1 <- p1 + geom_density(alpha=0.5)
    p1 <- p1 + geom_rug(aes(y=NULL, color=Novelty, position="jitter"))
    p1 <- p1 + scale_x_continuous(measure)

    p2 <- ggplot(data=molten, aes(x=Sample, y=value, group=Novelty, color=Novelty), y=value)
    p2 <- p2 + opts(title = paste(measure, name))
    p2 <- p2 + geom_smooth(alpha=0.5, aes(group=Novelty))
    p2 <- p2 + sampleTextLabel + sampleTextLabelScale
    p2 <- p2 + facet_grid(Novelty ~ ., scales="free")
    p2 <- p2 + xAxis
    
    distributePerSampleGraph(p1, p2)
  }
    
  # known / novel ratio by sample
  # TODO -- would ideally not conflate SNPs and Indels
  d = subset(metricsBySamples, Novelty == "all" & CompRod == "dbsnp")
  title <- opts(title = "Novelty rate by sample")

  # distribution
  p1 <- ggplot(data=d, aes(x=compRate))
  p1 <- p1 + geom_density(alpha=0.5)
  p1 <- p1 + geom_rug(aes(y=NULL, position="jitter"))
  p1 <- p1 + scale_x_continuous("Percent of variants in dbSNP")

  p2 <- ggplot(data=d, aes(x=Sample, y=compRate))
  p2 <- p2 + title
  p2 <- p2 + geom_smooth(alpha=0.5, aes(group=Novelty))
  p2 <- p2 + sampleTextLabel + sampleTextLabelScale
  p2 <- p2 + geom_rug(aes(x=NULL, position="jitter"))
  p2 <- p2 + xAxis
  p2 <- p2 + scale_y_continuous("Percent of variants in dbSNP")
  distributePerSampleGraph(p1, p2)

  for ( novelty in c("all", "known", "novel") ) {
    # TODO -- how can I color it as before?
    # TODO -- add marginal distributions?
    molten = melt(subset(metricsBySamples, Novelty==novelty), id.vars=c("Sample", "highlightTextSizes"), measure.vars=measures)
    p <- ggplot(data=molten, aes(x=Sample, y=value))
    p <- p + opts(title = paste(name, ":", novelty))
#    p <- p + scale_y_log10("Number of variants")
#    p <- p + geom_point(alpha=0.5, size=4)
    p <- p + sampleTextLabel + sampleTextLabelScale
    p <- p + facet_grid(variable ~ ., scales="free")
    # how do we remove the labels?
    p <- p + xAxis
    print(p)
  }
}

# -------------------------------------------------------
# Actually invoke the above plotting functions 
# -------------------------------------------------------

# load the data.
if ( onCMDLine || LOAD_DATA ) {
  metricsBySites <- createMetricsBySites(VariantEvalRoot)
  metricsBySamples <- createMetricsBySamples(VariantEvalRoot)
}

if ( ! is.na(outputPDF) ) {
  pdf(outputPDF, height=8.5, width=11)
} 

# Table of overall counts and quality
textplot(overallSummaryTable(metricsBySites), show.rownames=F)
title(paste("Summary metrics for project", ProjectName), cex=3)
# textplot(as.data.frame(sampleSummaryTable(metricsBySamples)), show.rownames=F)
# title(paste("Summary metrics per sample for project", ProjectName), cex=3)

summaryPlots(metricsBySites)
perSamplePlots(metricsBySamples)

if ( ! is.na(outputPDF) ) {
  dev.off()
} 

