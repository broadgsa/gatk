#!/bin/env Rscript

require("ggplot2")

args <- commandArgs(TRUE)
verbose = TRUE

inputDataFile = args[1]
onCmdLine = ! is.na(inputDataFile)

addEmpiricalPofG <- function(d) {
  r = c()
  #
  # TODO -- this is a really naive estimate of the accuracy, as it assumes the comp
  # track is perfect.  In reality the chip is at best Q30 accurate (replicate samples have
  # level than this level of concordance).  At low incoming confidence, we can effectively
  # ignore this term but when the incoming Q is near or above Q30 this approximation clearly
  # breaks down.
  #
  for ( i in 1:dim(d)[1] ) {
    row = d[i,]
    if ( row$pGGivenDType == "QofAAGivenD" ) v = row$HOM_REF
    if ( row$pGGivenDType == "QofABGivenD" ) v = row$HET
    if ( row$pGGivenDType == "QofBBGivenD" ) v = row$HOM_VAR
    r = c(r, v / row$Sum)
  }

  #print(length(r))
  d$EmpiricalPofG = r
  d$EmpiricalPofGQ = round(-10*log10(1-r))
  return(d)
}

genotypeCounts <- function(x) {
  type = unique(x$variable)[1]
  t = addmargins(table(x$comp))
  return(t)
}


digestTable <- function(inputDataFile) {
  d = subset(read.table(inputDataFile, header=T), rg != "ALL")
  d$technology <- factor(1, levels=c("HiSeq-paper", "GA2-1000G", "HiSeq-recent"))
  d$technology[grepl("ERR.*", d$rg)] <- "GA2-1000G"
  d$technology[grepl("20.*", d$rg)] <- "HiSeq-paper"
  d$technology[grepl("B00EG.*", d$rg)] <- "HiSeq-recent"
  print(summary(d$technology))
  
  eByComp = addEmpiricalPofG(ddply(d, .(rg, technology, pGGivenDType, pGGivenD), genotypeCounts))
  return(list(d=d, eByComp = eByComp))
  #countsByTech = addEmpiricalPofG(ddply(d, .(technology, pGGivenDType, pGGivenD), genotypeCounts))
}
  
writeMyTable <- function(t, name) {
  write.table(t,file=paste(inputDataFile, ".", name, ".tsv", sep=""))
}

if ( onCmdLine ) {
  r <- digestTable(inputDataFile)
  writeMyTable(r$eByComp, "eByComp")
}

