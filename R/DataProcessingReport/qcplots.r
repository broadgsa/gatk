.libPaths('/humgen/gsa-firehose2/pipeline/repositories/StingProduction/R/') 
.libPaths('~/Documents/Sting/R/') 

library(gsalib)
def.par <- par(no.readonly = TRUE)

titvplot<-function(current){
par(mfcol=c(1,2))
titvs<-c()
status<-c()
for(i in c(1:12)){ 
  load(sprintf("%sexome.%i", path, i)); 
  info<-subset(data$TiTvVariantEvaluator, Sample!="all")
  titvs<-c(titvs, info$tiTvRatio)
  status<-c(status, info$Novelty)
  print(length(titvs))
  print(length(status))
  }
print(length(unique(current$TiTvVariantEvaluator$Sample))-1)

length(unique(current$TiTvVariantEvaluator$Sample))-1+length(titvs[which(status=="novel")])->nvalues
print(length(titvs[which(status=="novel")]))
print(nvalues)
plot(current$TiTvVariantEvaluator$tiTvRatio[which(current$TiTvVariantEvaluator$Sample!="all" & current$TiTvVariantEvaluator$Novelty=="novel")], xlim=c(0,nvalues), ylim=c(0,4), col="red", main="Current samples compared to previous samples from 12 sets", ylab="Per sample Ti/Tv", xlab="sample")
points(current$TiTvVariantEvaluator$tiTvRatio[which(current$TiTvVariantEvaluator$Sample!="all" & current$TiTvVariantEvaluator$Novelty=="known")], col="blue")
points(current$TiTvVariantEvaluator$tiTvRatio[which(current$TiTvVariantEvaluator$Sample!="all" & current$TiTvVariantEvaluator$Novelty=="all")],  col="black")
points(c(length(unique(current$TiTvVariantEvaluator$Sample)):nvalues), titvs[which(status=="novel")], pch=16, col="red")
points(c(length(unique(current$TiTvVariantEvaluator$Sample)):nvalues), titvs[which(status=="known")], pch=16, col="blue")
points(c(length(unique(current$TiTvVariantEvaluator$Sample)):nvalues), titvs[which(status=="all")], pch=16, col="black")

legend("bottomleft", col=c("red", "blue", "black"), pch=c(1,1,1,16,16, 16),legend=c("novel variants:current set", "known variants:current set", "all varaints:current set", "novel variants:previous sets", "known variants:previous sets", "all variants: previous sets"))  
weirdos<-which(current$TiTvVariantEvaluator$Sample %in% current$TiTvVariantEvaluator$Sample[which(current$TiTvVariantEvaluator$tiTvRatio <2.0)])
if(length(weirdos)>0){
  text(weirdos[c(1:(length(weirdos)/3))],current$TiTvVariantEvaluator$tiTvRatio[weirdos], labels=current$TiTvVariantEvaluator$Sample[weirdos], pos=4, cex=.7, col="hot pink")
}

boxplot(current$TiTvVariantEvaluator$tiTvRatio[which(current$TiTvVariantEvaluator$Sample!="all" & current$TiTvVariantEvaluator$Novelty=="novel")],titvs[which(status=="novel")], current$TiTvVariantEvaluator$tiTvRatio[which(current$TiTvVariantEvaluator$Sample!="all" & current$TiTvVariantEvaluator$Novelty=="known")],titvs[which(status=="known")], current$TiTvVariantEvaluator$tiTvRatio[which(current$TiTvVariantEvaluator$Sample!="all" & current$TiTvVariantEvaluator$Novelty=="all")], titvs[which(status=="all")], col=rep(c("red", "blue", "black"), each=2), main="Current v. Previous per sample Ti/TV", xlab="Sample Sets",ylab="Ti/Tv per sample", xaxt="n" )
axis(side=1, at=c(1:6)-.2, labels=rep(c("current", "previous"), 3), cex.axis=.7)
legend("bottomleft",legend=c("novel", "known", "all"), fill=c("red", "blue", "black"))
if(length(weirdos)>0){
text(rep(c(5,3,1), each=(length(weirdos)/3)),current$TiTvVariantEvaluator$tiTvRatio[weirdos], labels=current$TiTvVariantEvaluator$Sample[weirdos], pos=4, cex=.7, col="hot pink")
}
par(def.par)#- reset to default

}



variantplots<-function(current){
  par(mfcol=c(1,2))

variants<-c()
status<-c()
for(i in c(1:12)){ 
  load(sprintf("%s/exome.%i", path, i)); 
  info<-subset(data$CountVariants, Sample!="all")
  variants<-c(variants, info$nSNPs)
  status<-c(status, info$Novelty)
  }

length(unique(current$CountVariants$Sample))-1+length(variants[which(status=="novel")])->nvalues
plot(current$CountVariants$nSNPs[which(current$CountVariants$Sample!="all" & current$CountVariants$Novelty=="novel")], xlim=c(0,nvalues), ylim=c(1,25000), log="y", col="red", main="Current samples compared to previous samples from 12 sets", ylab="Per sample #SNPs", xlab="sample")
points(current$CountVariants$nSNPs[which(current$CountVariants$Sample!="all" & current$CountVariants$Novelty=="known")], col="blue")
points(current$CountVariants$nSNPs[which(current$CountVariants$Sample!="all" & current$CountVariants$Novelty=="all")],  col="black")
points(c(length(unique(current$CountVariants$Sample)):nvalues), variants[which(status=="novel")], pch=16, col="red")
points(c(length(unique(current$CountVariants$Sample)):nvalues), variants[which(status=="known")], pch=16, col="blue")
points(c(length(unique(current$CountVariants$Sample)):nvalues), variants[which(status=="all")], pch=16, col="black")

legend("bottomleft", col=c("red", "blue", "black"), pch=c(1,1,1,16,16, 16),legend=c("novel variants:current set", "known variants:current set", "all varaints:current set", "novel variants:previous sets", "known variants:previous sets", "all variants: previous sets"))  

weirdos<-which(current$CountVariants$Sample %in% current$TiTvVariantEvaluator$Sample[which(current$TiTvVariantEvaluator$tiTvRatio <2.0)])
if(length(weirdos)>0){

text(weirdos[c(1:(length(weirdos)/3))],current$CountVariants$nSNPs[weirdos], labels=current$CountVariants$Sample[weirdos], pos=4, cex=.7, col="hot pink")
}

boxplot(current$CountVariants$nSNPs[which(current$CountVariants$Sample!="all" & current$CountVariants$Novelty=="novel")],variants[which(status=="novel")], current$CountVariants$nSNPs[which(current$CountVariants$Sample!="all" & current$CountVariants$Novelty=="known")],variants[which(status=="known")], current$CountVariants$nSNPs[which(current$CountVariants$Sample!="all" & current$CountVariants$Novelty=="all")], variants[which(status=="all")], col=rep(c("red", "blue", "black"), each=2), main="Current v. Previous per sample #SNPs", xlab="Sample Sets",ylab="SNPs per sample", xaxt="n", ylim=c(10,25000), log="y")
axis(side=1, at=c(1:6)-.2, labels=rep(c("current", "previous"), 3), cex.axis=.7)
 if(length(weirdos)>0){

 text(rep(c(5,3,1), each=(length(weirdos)/3)),current$CountVariants$nSNPs[weirdos], labels=current$CountVariants$Sample[weirdos], pos=4, cex=.7, col="hot pink")
}
legend("topleft",legend=c("novel", "known", "all"), fill=c("red", "blue", "black"))
par(def.par)#- reset to default

}

heteroplots<-function(current){
  par(mfcol=c(1,2))

hets<-c()
status<-c()
for(i in c(1:12)){ 
  load(sprintf("%s/exome.%i", path, i)); 
  info<-subset(data$CountVariants, Sample!="all")
  hets<-c(hets, info$heterozygosity)
  status<-c(status, info$Novelty)
  }

length(unique(current$CountVariants$Sample))-1+length(hets[which(status=="novel")])->nvalues
plot(current$CountVariants$heterozygosity[which(current$CountVariants$Sample!="all" & current$CountVariants$Novelty=="novel")], xlim=c(0,nvalues), ylim=c(-0.0005, 0.0005), col="red", main="Current samples compared to previous samples from 12 sets", ylab="Per sample heterozygosity", xlab="sample")
points(current$CountVariants$heterozygosity[which(current$CountVariants$Sample!="all" & current$CountVariants$Novelty=="known")], col="blue")
points(current$CountVariants$heterozygosity[which(current$CountVariants$Sample!="all" & current$CountVariants$Novelty=="all")],  col="black")
points(c(length(unique(current$CountVariants$Sample)):nvalues), hets[which(status=="novel")], pch=16, col="red")
points(c(length(unique(current$CountVariants$Sample)):nvalues), hets[which(status=="known")], pch=16, col="blue")
points(c(length(unique(current$CountVariants$Sample)):nvalues), hets[which(status=="all")], pch=16, col="black")

legend("bottomleft", col=c("red", "blue", "black"), pch=c(1,1,1,16,16, 16),legend=c("novel variants:current set", "known variants:current set", "all varaints:current set", "novel variants:previous sets", "known variants:previous sets", "all variants: previous sets"))  

weirdos<-which(current$CountVariants$Sample %in% current$TiTvVariantEvaluator$Sample[which(current$TiTvVariantEvaluator$tiTvRatio <2.0)])
 if(length(weirdos)>0){
text(weirdos[c(1:(length(weirdos)/3))],current$CountVariants$heterozygosity[weirdos], labels=current$CountVariants$Sample[weirdos], pos=4, cex=.7, col="hot pink")
}

boxplot(current$CountVariants$heterozygosity[which(current$CountVariants$Sample!="all" & current$CountVariants$Novelty=="novel")],hets[which(status=="novel")], current$CountVariants$heterozygosity[which(current$CountVariants$Sample!="all" & current$CountVariants$Novelty=="known")],hets[which(status=="known")], current$CountVariants$heterozygosity[which(current$CountVariants$Sample!="all" & current$CountVariants$Novelty=="all")], hets[which(status=="all")], col=rep(c("red", "blue", "black"), each=2), main="Current v. Previous per sample #Heterozygousity", xlab="Sample Sets",ylab="Heterozygousity per sample", xaxt="n")
axis(side=1, at=c(1:6)-.2, labels=rep(c("current", "previous"), 3), cex.axis=.7)
if(length(weirdos)>0){

text(rep(c(5,3,1), each=(length(weirdos)/3)),current$CountVariants$heterozygosity[weirdos], labels=current$CountVariants$Sample[weirdos], pos=4, cex=.7, col="hot pink")
}
legend("topleft",legend=c("novel", "known", "all"), fill=c("red", "blue", "black"))
par(def.par)#- reset to default

}

novelAC<-function(current){
ACs<-sort(current$SimpleMetricsByAC.metrics$AC[which(current$SimpleMetricsByAC.metrics$Novelty=="novel")])
orderbyAC<-order(current$SimpleMetricsByAC.metrics$AC[which(current$SimpleMetricsByAC.metrics$Novelty=="novel")])
varbyAC<-current$SimpleMetricsByAC.metrics$n[which(current$SimpleMetricsByAC.metrics$Novelty=="novel")][orderbyAC]
plot(ACs, varbyAC, type="l", log="xy", lwd=4, col="dark red", main="Novel AC", ylab="# variants (log scale)", xlab="AC (log scale)")

for(i in c(1:12)){ 
  load(sprintf("%s/exome.%i", path, i)); 
  info<-data$SimpleMetricsByAC.metrics
  ACs<-sort(info$AC[which(info$Novelty=="novel")])
  orderbyAC<-order(info$AC[which(info$Novelty=="novel")])
  varbyAC<-info$n[which(info$Novelty=="novel")][orderbyAC]

  lines(ACs, varbyAC, col="red")
}

legend("topright",legend=c("current", "previous"), lwd=c(4,1), col=c("dark red", "red"))
}

knownAC<-function(current){
ACs<-sort(current$SimpleMetricsByAC.metrics$AC[which(current$SimpleMetricsByAC.metrics$Novelty=="known")])
orderbyAC<-order(current$SimpleMetricsByAC.metrics$AC[which(current$SimpleMetricsByAC.metrics$Novelty=="known")])
varbyAC<-current$SimpleMetricsByAC.metrics$n[which(current$SimpleMetricsByAC.metrics$Novelty=="known")][orderbyAC]
plot(ACs, varbyAC, type="l", log="xy", lwd=4, col="dark blue", main="Known AC", ylab="# variants (log scale)", xlab="AC (log scale)")

for(i in c(1:12)){ 
  load(sprintf("%s/exome.%i", path, i)); 
  info<-data$SimpleMetricsByAC.metrics
  ACs<-sort(info$AC[which(info$Novelty=="known")])
  orderbyAC<-order(info$AC[which(info$Novelty=="known")])
  varbyAC<-info$n[which(info$Novelty=="known")][orderbyAC]
  lines(ACs, varbyAC, col="light blue")
}

legend("topright",legend=c("current", "previous"), lwd=c(4,1), col=c("dark blue", "light blue"))
}

AllAC<-function(current){
ACs<-sort(current$SimpleMetricsByAC.metrics$AC[which(current$SimpleMetricsByAC.metrics$Novelty=="all")])
orderbyAC<-order(current$SimpleMetricsByAC.metrics$AC[which(current$SimpleMetricsByAC.metrics$Novelty=="all")])
varbyAC<-current$SimpleMetricsByAC.metrics$n[which(current$SimpleMetricsByAC.metrics$Novelty=="all")][orderbyAC]
plot(ACs, varbyAC, type="l", log="xy", lwd=4, col="Black", main="All AC", ylab="# variants (log scale)", xlab="AC (log scale)")

for(i in c(1:12)){ 
  load(sprintf("%s/exome.%i", path, i)); 
  info<-data$SimpleMetricsByAC.metrics
  ACs<-sort(info$AC[which(info$Novelty=="all")])
  orderbyAC<-order(info$AC[which(info$Novelty=="all")])
  varbyAC<-info$n[which(info$Novelty=="all")][orderbyAC]

  lines(ACs, varbyAC, col="dark grey")
}

legend("topright",legend=c("current", "previous"), lwd=c(4,1), col=c("black", "dark grey"))
}

