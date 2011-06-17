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
  preQCFile <- "~/Desktop/broadLocal/GATK/trunk/qcTestData/GoT2D_exomes_batch_005_per_sample_metrics.tsv"
  VariantEvalRoot <- "~/Desktop/broadLocal/GATK/trunk/qcTestData/GoT2D_exomes_batch_005.cleaned.snps_and_indels.filtered.annotated"
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

createMetricsBySites <- function(VariantEvalRoot, PreQCMetrics) {
  # Metrics by sites:
  #  bySite -> counts of SNPs and Indels by novelty, with expectations
  #  byAF -> snps and indels (known / novel)
  r = list( bySite = expandVEReport(gsa.read.gatkreport(paste(VariantEvalRoot, ".summary.eval", sep=""))),
               byAF = gsa.read.gatkreport(paste(VariantEvalRoot, ".byAF.eval", sep="")))
  r$byAF$CountVariants$nIndels = r$byAF$CountVariants$nInsertions + r$byAF$CountVariants$nDeletions

  nChrom = 1 / min(r$byAF$TiTvVariantEvaluator$AlleleFrequency[r$byAF$TiTvVariantEvaluator$AlleleFrequency > 0])
  r$byAF$TiTvVariantEvaluator$nSNPs = r$byAF$TiTvVariantEvaluator$nTi + r$byAF$TiTvVariantEvaluator$nTv
  r$byAF$CountVariants$AC = r$byAF$CountVariants$AlleleFrequency * nChrom
  r$byAF$TiTvVariantEvaluator$AC = r$byAF$TiTvVariantEvaluator$AlleleFrequency * nChrom
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

summaryPlots <- function(metricsBySites) {
  name = "SNP and Indel count by novelty and allele frequency" 
  molten = melt(subset(metricsBySites$byAF$CountVariants, FunctionalClass == "all" & Novelty != "all" & AC > 0), id.vars=c("Novelty", "AC"), measure.vars=c(c("nSNPs", "nIndels")))
  p <- ggplot(data=molten, aes(x=AC, y=value+1, color=Novelty, fill=Novelty), group=variable)
  p <- p + opts(title = name)
  p <- p + scale_y_log10("Number of variants")
  p <- p + geom_point(alpha=0.5, size=3)
  p <- p + geom_line(size=1)
  p <- p + facet_grid(variable ~ ., scales="free")
  p <- p + scale_x_continuous("Allele count (AC)")
  print(p)
  p <- p + scale_x_log10("Allele count (AC)")
  print(p)

  # Counts vs. Allele frequency 
  name = "Variant counts by allele count"
  for ( measure in c("nSNPs", "nIndels")) {
    molten = melt(subset(metricsBySites$byAF$CountVariants, FunctionalClass == "all" & AC > 0), id.vars=c("Novelty", "AC"), measure.vars=c(measure))
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
  byAFNoAll = subset(metricsBySites$byAF$TiTvVariantEvaluator, Novelty != "all"  & FunctionalClass == "all" & AC > 0)
  p <- ggplot(data=byAFNoAll, aes(x=AC, y=tiTvRatio, color=Novelty))
  p <- p + scale_y_continuous("Transition / transversion ratio", limits=c(0,4))
  p <- p + opts(title = name)
  p <- p + geom_smooth(size=2)
  p <- p + geom_point(aes(size=log10(nSNPs), weight=nSNPs), alpha=0.5)
  p <- p + scale_x_continuous("Allele count (AC)")
  print(p)
  p <- p + scale_x_log10("Allele count (AC)")
  print(p)
  
  # SNPs to indels ratio by allele frequency
  name = "SNPs to indels ratio by allele frequency" 
  metricsBySites$byAF$CountVariants$SNP.Indel.Ratio = metricsBySites$byAF$CountVariants$nSNPs / metricsBySites$byAF$CountVariants$nIndels
  metricsBySites$byAF$CountVariants$SNP.Indel.Ratio[metricsBySites$byAF$CountVariants$nIndels == 0] = NaN
  p <- ggplot(data=subset(metricsBySites$byAF$CountVariants, FunctionalClass == "all" & Novelty == "all"), aes(x=AC, y=SNP.Indel.Ratio))
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

  #r = merge(merge(preQCMetrics, byAFEval$TiTvVariantEvaluator), byAFEval$CountVariants)
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
  myRug <- geom_rug(position="jitter")
  #myRug <- geom_rug(aes(x=NULL), position="jitter")

  measures = c("nSNPs", "tiTvRatio", "nSingletons", "nIndels", "nInsertions", "nDeletions", "deletionInsertionRatio")
  name = "by sample"
  for ( measure in measures ) {
    molten = melt(metricsBySamples, id.vars=c("Novelty", "Sample", "highlightTextSizes"), measure.vars=c(measure))

    # distribution
    p <- ggplot(data=molten, aes(x=value, group=Novelty, fill=Novelty))
    p <- p + opts(title = paste(measure, name))
    p <- p + geom_density(alpha=0.5)
    p <- p + geom_rug(aes(y=NULL, color=Novelty, position="jitter"))
    p <- p + scale_x_continuous(measure)
    # how do we remove the labels?
    print(p)

    p <- ggplot(data=molten, aes(x=Sample, y=value, group=Novelty, color=Novelty), y=value)
    p <- p + opts(title = paste(measure, name))
    p <- p + geom_smooth(alpha=0.5, aes(group=Novelty))
    p <- p + sampleTextLabel + sampleTextLabelScale
    p <- p + myRug
    p <- p + facet_grid(Novelty ~ ., scales="free")
    # how do we remove the labels?
    p <- p + xAxis
    print(p)
  }

  # known / novel ratio by sample
  # TODO -- would ideally not conflate SNPs and Indels
  d = subset(metricsBySamples, Novelty == "all" & CompRod == "dbsnp")
  title <- opts(title = "Novelty rate by sample")

  # distribution
  p <- ggplot(data=d, aes(x=compRate))
  p <- p + title
  p <- p + geom_density(alpha=0.5)
  p <- p + geom_rug(aes(y=NULL, position="jitter"))
  p <- p + scale_x_continuous("Percent of variants in dbSNP")
  # how do we remove the labels?
  print(p)

  p <- ggplot(data=d, aes(x=Sample, y=compRate))
  p <- p + title
  p <- p + geom_smooth(alpha=0.5, aes(group=Novelty))
  p <- p + sampleTextLabel + sampleTextLabelScale
  p <- p + geom_rug(aes(x=NULL, position="jitter"))
  #p <- p + myRug
  # how do we remove the labels?
  p <- p + xAxis
  print(p)

  for ( novelty in c("all", "known", "novel") ) {
    # TODO -- how can I color it as before?
    # TODO -- add marginal distributions?
    molten = melt(subset(metricsBySamples, Novelty==novelty), id.vars=c("Sample", "highlightTextSizes"), measure.vars=measures)
    p <- ggplot(data=molten, aes(x=Sample, y=value))
    p <- p + opts(title = paste(name, ":", novelty))
#    p <- p + scale_y_log10("Number of variants")
#    p <- p + geom_point(alpha=0.5, size=4)
    p <- p + sampleTextLabel + sampleTextLabelScale
    p <- p + myRug
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
textplot(as.data.frame(summaryTable(metricsBySites)), show.rownames=F)
title(paste("Overall summary metrics for project", ProjectName), cex=3)
 
summaryPlots(metricsBySites)
perSamplePlots(metricsBySamples)

if ( ! is.na(outputPDF) ) {
  dev.off()
} 


