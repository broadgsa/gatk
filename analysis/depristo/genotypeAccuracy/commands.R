require("lattice")
require("ggplot2")

READ_DATA = F

if ( READ_DATA ) {
  d = read.table("~/Dropbox/Analysis/genotypeAccuracy/cgl.hiseq.table", header=T)
  #d = read.table("~/Desktop/broadLocal/GATK/trunk/foo", header=T)
}

moltenCD = d
#moltenCD = melt(d, id.vars=c("comp", "rg"), measure.vars=c("QofAAGivenD", "QofABGivenD", "QofBBGivenD"))
#moltenCD$log10value = round(-10*log10(1-10^moltenCD$value))
  
genotypeCounts <- function(x) {
  #print(table(x$comp))
  type = unique(x$variable)[1]
  #print(type)
  t = addmargins(table(x$comp))
  return(t)
}

addEmpiricalPofG <- function(d) {
  r = c()
  for ( i in 1:dim(d)[1] ) {
    row = d[i,]
    #print(row)
    if ( row$pGGivenDType == "QofAAGivenD" ) v = row$HOM_REF
    if ( row$pGGivenDType == "QofABGivenD" ) v = row$HET
    if ( row$pGGivenDType == "QofBBGivenD" ) v = row$HOM_VAR
    #print(v)
    #print(row$Sum)
    r = c(r, v / row$Sum)
    #print(r)
  }

  #print(length(r))
  d$EmpiricalPofG = r
  d$EmpiricalPofGQ = round(-10*log10(1-r))
  return(d)
}

eByComp <- addEmpiricalPofG(ddply(moltenCD, .(rg, pGGivenDType, pGGivenD), genotypeCounts))
print(subset(eByComp, EmpiricalPofGQ < Inf))

goodEByComp = subset(eByComp, Sum > 10 & EmpiricalPofGQ < Inf)

print(qplot(pGGivenD, EmpiricalPofGQ, data=subset(goodEByComp, rg != "ALL"), size=log10(Sum), color=pGGivenDType, geom=c("point", "smooth"), group=pGGivenDType, xlim=c(0,40), ylim=c(0,40)) + geom_abline(slope=1, linetype=2))
print(qplot(pGGivenD, EmpiricalPofGQ, data=subset(goodEByComp, rg != "ALL"), facets = . ~ pGGivenDType, color=rg, geom=c("line"), group=rg, xlim=c(0,40), ylim=c(0,40)) + geom_abline(slope=1, linetype=2))



