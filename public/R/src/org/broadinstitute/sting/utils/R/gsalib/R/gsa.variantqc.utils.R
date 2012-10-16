library(gplots)
library(ggplot2)
library(tools)

# -------------------------------------------------------
# Utilities for displaying multiple plots per page
# -------------------------------------------------------

distributeGraphRows <- function(graphs, heights = c()) {
  # Viewport layout 2 graphs top to bottom with given relative heights
  #
  #
  if (length(heights) == 0) {
    heights <- rep.int(1, length(graphs))
  }
  heights <- heights[!is.na(graphs)]
  graphs <- graphs[!is.na(graphs)]
  numGraphs <- length(graphs)
  Layout <- grid.layout(nrow = numGraphs, ncol = 1, heights=heights)
  grid.newpage()
  pushViewport(viewport(layout = Layout))
  subplot <- function(x) viewport(layout.pos.row = x, layout.pos.col = 1)
  for (i in 1:numGraphs) {
    print(graphs[[i]], vp = subplot(i))
  }
}

distributeLogGraph <- function(graph, xName) {
  continuousGraph <- graph + scale_x_continuous(xName)
  logGraph <- graph + scale_x_log10(xName) + opts(title="")
  distributeGraphRows(list(continuousGraph, logGraph))
}

distributePerSampleGraph <- function(perSampleGraph, distGraph, ratio=c(2,1)) {
  distributeGraphRows(list(perSampleGraph, distGraph), ratio)
}

removeExtraStrats <- function(variantEvalDataFrame, moreToRemove=c()) {
  # Remove the standard extra stratification columns FunctionalClass, Novelty, and others in moreToRemove from the variantEvalDataFrame
  #
  # Only keeps the column marked with "all" for each removed column
  #
  for ( toRemove in c("FunctionalClass", "Novelty", moreToRemove) ) {
    if (toRemove %in% colnames(variantEvalDataFrame)) {
      variantEvalDataFrame <- variantEvalDataFrame[variantEvalDataFrame[[toRemove]] == "all",]
    }
  }
  variantEvalDataFrame    
}

openPDF <- function(outputPDF) {
  # Open the outputPDF file with standard dimensions, if outputPDF is not NA
  if ( ! is.na(outputPDF) ) {
    pdf(outputPDF, height=8.5, width=11)
  }
}

closePDF <- function(outputPDF) {
  # close the outputPDF file if not NA, and try to compact the PDF if possible
  if ( ! is.na(outputPDF) ) {
    dev.off()
    if (exists("compactPDF")) {
      print("compacting PDF")
      compactPDF(outputPDF)
    }
  }
}

makeRatioDataFrame <- function(ACs, num, denom, widths = NULL) {
  if ( is.null(widths) ) widths <- rep(1, length(ACs))
  
  value = NULL
  titv <- data.frame(AC=ACs, width = widths, num=num, denom = denom, ratio = num / denom)
}

.reduceACs <- function(binWidthForAC, ACs) {
  # computes data structures necessary to reduce the full range of ACs
  #
  # binWidthForAC returns the number of upcoming bins that should be merged into 
  # that AC bin.  ACs is a vector of all AC values from 0 to 2N that should be 
  # merged together
  #
  # Returns a list containing the reduced ACs starts, their corresponding widths,
  # and a map from original ACs to their new ones (1 -> 1, 2 -> 2, 3 -> 2, etc)
  maxAC <- max(ACs)
  newACs <- c()
  widths <- c()
  newACMap <- c()
  ac <- 0
  while ( ac < maxAC ) {
    newACs <- c(newACs, ac)
    width <- binWidthForAC(ac)
    widths <- c(widths, width)
    newACMap <- c(newACMap, rep(ac, width))
    ac <- ac + width
  }
  list(ACs = newACs, widths=widths, newACMap = newACMap)
}

# geometricACs <- function(k, ACs) {
#   nBins <- round(k * log10(max(ACs)))
#   
#   binWidthForAC <- function(ac) {
#     max(ceiling(ac / nBins), 1)
#   }
#   
#   return(reduceACs(binWidthForAC, ACs))
# }

reduce.AC.on.LogLinear.intervals <- function(scaleFactor, ACs) {
  # map the full range of AC values onto a log linear scale
  #
  # Reduce the full AC range onto one where the width of each new AC increases at a rate of
  # 10^scaleFactor in size with growing AC values.  This is primarily useful for accurately
  # computing ratios or other quantities by AC that aren't well determined when the AC 
  # values are very large
  #
  # Returns a list containing the reduced ACs starts, their corresponding widths,
  # and a map from original ACs to their new ones (1 -> 1, 2 -> 2, 3 -> 2, etc)
  maxAC <- max(ACs)
  afs <- ACs / maxAC
  breaks <- 10^(seq(-4, -1, scaleFactor))
  widths <- c()
  lastBreak <- 1
  for ( i in length(breaks):1 ) {
    b <- breaks[i]
    width <- sum(afs < lastBreak & afs >= b)
    widths <- c(widths, width)
    lastBreak <- b
  }
  widths <- rev(widths)
  
  binWidthForAC <- function(ac) {
    af <- ac / maxAC
    value = 1
    for ( i in length(breaks):1 )
      if ( af >= breaks[i] ) {
        value = widths[i]
        break
      }
    
    return(value)
  }
  
  return(.reduceACs(binWidthForAC, ACs))
}

.remapACs <- function(remapper, k, df) {
  newACs <- remapper(k, df$AC)
  
  n = length(newACs$ACs)
  num = rep(0, n)
  denom = rep(0, n)
  for ( i in 1:dim(df)[1] ) {
    rowI = df$AC == i
    row = df[rowI,]
    newAC = newACs$newACMap[row$AC]
    newRowI = newACs$ACs == newAC
    num[newRowI] = num[newRowI] + df$num[rowI]
    denom[newRowI] = denom[newRowI] + df$denom[rowI]
  }
  
  newdf <- makeRatioDataFrame(newACs$ACs, num, denom, newACs$widths )
  newdf
}

compute.ratio.on.LogLinear.AC.intervals <- function(ACs, num, denom, scaleFactor = 0.1) {
  df = makeRatioDataFrame(ACs, num, denom, 1)
  return(.remapACs(reduce.AC.on.LogLinear.intervals, scaleFactor, df))
}

plotVariantQC <- function(metrics, measures, requestedStrat = "Sample", 
                          fixHistogramX=F, anotherStrat = NULL, nObsField = "n_indels", 
                          onSamePage=F, facetVariableOnXPerSample = F, facetVariableOnXForDist = T, 
                          moreTitle="", note = NULL) {
  metrics$strat = metrics[[requestedStrat]]
  
  otherFacet = "."
  id.vars = c("strat", "nobs")
  metrics$nobs <- metrics[[nObsField]]
  
  # keep track of the other strat and it's implied facet value
  if (! is.null(anotherStrat)) { 
    id.vars = c(id.vars, anotherStrat)
    otherFacet = anotherStrat
  }
  
  molten <- melt(metrics, id.vars=id.vars, measure.vars=c(measures))
  perSampleGraph <- ggplot(data=molten, aes(x=strat, y=value, group=variable, color=variable, fill=variable))

  # create the title
  titleText=paste(paste(paste(measures, collapse=", "), "by", requestedStrat), moreTitle)
  if ( !is.null(note) ) {
    titleText=paste(titleText, note, sep="\n")
  }
  paste(titleText)
  title <- opts(title=titleText)
  
  determineFacet <- function(onX) {
    if ( onX ) { 
      paste(otherFacet, "~ variable")
    } else {
      paste("variable ~", otherFacet)
    }
  }
  
  sampleFacet = determineFacet(facetVariableOnXPerSample)
  distFacet   = determineFacet(facetVariableOnXForDist)
  
  if ( requestedStrat == "Sample" ) {
    perSampleGraph <- perSampleGraph + geom_text(aes(label=strat), size=1.5) + geom_blank() # don't display a scale
    perSampleGraph <- perSampleGraph + scale_x_discrete("Sample (ordered by nSNPs)")
  } else { # by AlleleCount
    perSampleGraph <- perSampleGraph + geom_point(aes(size=log10(nobs))) #+ geom_smooth(aes(weight=log10(nobs)))
    perSampleGraph <- perSampleGraph + scale_x_log10("AlleleCount")
  }    
  perSampleGraph <- perSampleGraph + ylab("Variable value") + title
  perSampleGraph <- perSampleGraph + facet_grid(sampleFacet, scales="free")
  
  nValues = length(unique(molten$value))
  if (nValues > 2) {
    if ( requestedStrat == "Sample" ) {
      distGraph <- ggplot(data=molten, aes(x=value, group=variable, fill=variable))
    } else {
      distGraph <- ggplot(data=molten, aes(x=value, group=variable, fill=variable, weight=nobs))
    }
    distGraph <- distGraph + geom_histogram(aes(y=..ndensity..))
    distGraph <- distGraph + geom_density(alpha=0.5, aes(y=..scaled..))
    distGraph <- distGraph + geom_rug(aes(y=NULL, color=variable, position="jitter"))
    scale = "free"
    if ( fixHistogramX ) scale = "fixed"
    distGraph <- distGraph + facet_grid(distFacet, scales=scale)
    distGraph <- distGraph + ylab("Relative frequency")
    distGraph <- distGraph + xlab("Variable value (see facet for variable by color)")
    distGraph <- distGraph + opts(axis.text.x=theme_text(angle=-45)) # , legend.position="none")
  } else {
    distGraph <- NA
  }
  
  if ( onSamePage ) {
    suppressMessages(distributePerSampleGraph(perSampleGraph, distGraph))
  } else {
    suppressMessages(print(perSampleGraph))
    suppressMessages(print(distGraph + title))
  }
}
