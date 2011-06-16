require("lattice")
require("ggplot2")
require("splines")

ymax = xmax = 30
HAVE_RAW_DATA = F
if ( HAVE_RAW_DATA ) {
  inputDataFile = "~/Dropbox/Analysis/genotypeAccuracy/NA12878.hm3.vcf.cgl.table"
  #inputDataFile = "~/Dropbox/Analysis/genotypeAccuracy/cgl.table.gz"
  r <- digestTable(inputDataFile)
  d = r$d
  eByComp = r$eByComp
  countsByTech = addEmpiricalPofG(ddply(d, .(ref, alt, technology, pGGivenDType, pGGivenD), genotypeCounts))
  print(qplot(pGGivenD, EmpiricalPofGQ, data=subset(countsByTech, technology=="HiSeq-paper" & pGGivenDType == "QofABGivenD"), facets = alt ~ ref, color=alt, geom=c("point"), group=alt, xlim=c(0,xmax), ylim=c(0,ymax))
   + geom_abline(slope=1, linetype=2))
   #   + geom_smooth(se=T, size=1.5, aes(weight=Sum)))
} else {
  eByComp = read.table("~/Dropbox/GSA members/Analysis/genotypeAccuracy/NA12878.hm3.vcf.cgl.table.eByComp.tsv", header=T)
}

#print(subset(countsByTech, pGGivenD > 18 & pGGivenD < 22 & pGGivenDType == "QofABGivenD"))
#print(subset(eByComp, EmpiricalPofGQ < Inf))

goodEByComp = subset(eByComp, Sum > 10 & EmpiricalPofGQ < Inf)

print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, size=log10(Sum), facets = pGGivenDType ~ technology, color=pGGivenDType, geom=c("point", "smooth"), group=pGGivenDType, xlim=c(0,xmax), ylim=c(0,ymax)) + geom_abline(slope=1, linetype=2))

print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ technology, color=rg, geom=c("blank"), group=rg, xlim=c(0,xmax), ylim=c(0,ymax)) 
  + geom_abline(slope=1, linetype=2)
  + geom_smooth(se=F, aes(weight=Sum)))

print(qplot(pGGivenD, pGGivenD - EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ technology, color=rg, geom=c("blank"), group=rg, xlim=c(0,xmax), ylim=c(-10,10)) 
  + geom_abline(slope=0, linetype=2)
  + geom_smooth(se=F, method=lm, formula = y ~ ns(x,1), aes(weight=Sum)))

# By tech
print(qplot(pGGivenD, EmpiricalPofGQ, data=goodEByComp, facets = pGGivenDType ~ ., color=technology, geom=c("blank"), group=technology, xlim=c(0,xmax), ylim=c(0,ymax)) 
+ geom_abline(slope=1, linetype=2)
+ geom_smooth(se=T, size=1.5, aes(weight=Sum)))

