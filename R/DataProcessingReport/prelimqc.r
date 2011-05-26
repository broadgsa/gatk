#preqc.r
library(gplots)
.libPaths('/humgen/gsa-pipeline/.repository/R/') 
library(gsalib)

cmdargs = gsa.getargs(
    list(
        tsv = list(value=NA, doc="pipeline tsv file"),
        qcout=list(value=NA, doc="path to output root")
    ),
    doc="Creates a tearsheet"
);

read.delim(cmdargs$tsv, header=FALSE)->settable

 squids<-unique(settable[,1])
print(paste(nrow(settable), "samples in tsv"))
lane<-data.frame()
samp<-data.frame()
for(squid in squids){
  gsa.read.squidmetrics(squid, TRUE)->lanemetrics 
  print(paste("Got lane metrics for", squid))
  addlanes<-lanemetrics[which(lanemetrics$"External ID" %in% settable[,2]),]
  gsa.read.squidmetrics(squid, FALSE)->samplemetrics 
  print(paste("Got sample metrics for", squid))
  addsamps<-samplemetrics[which(samplemetrics$Sample %in% settable[,2]),]
  lane<-rbind(lane, addlanes)
  samp<-rbind(samp, addsamps)
}

print(paste(nrow(samp), "samples in samp"))
print(paste(length(unique(lane$"External ID")), "samples in lane"))

print(paste(setdiff(settable[,2], samp$Sample), "do not overlap between samp and tsv"))
print(paste(setdiff(settable[,2], lane$"External ID"), "do not overlap between lane and tsv"))
print(paste(setdiff(samp$Sample, lane$"External ID"), "do not overlap between lane and samp"))

missingSamp<-setdiff(settable[,2], samp$Sample) 
missingLane<-setdiff(settable[,2], lane$"External ID")

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

  compsamp=d2[which(d2$"Bait Set" %in% samp$"Bait Set"),]
  complane=d[which(d$"Bait Set" %in% lane$"Bait Set"),]



pdf(paste(cmdargs$qcout, "pdf", sep="."), width=11, height=8.5)

plot(samp$"Target Bases 20x %", main="Coverage to 20x", ylab="% Targets Covered to 20x", xlab="Sample", ylim=c(0,100))
abline(h=80, lty=2)
legend("bottomright", lty=2, legend="80% coverage to 20x")
lowcoverage<-samp$Sample[which(samp$"Target Bases 20x %"<80)]
if(length(lowcoverage)>0){
text(which(samp$"Target Bases 20x %"<80),samp$"Target Bases 20x %"[which(samp$"Target Bases 20x %"<80)], labels=samp$Sample[which(samp$"Target Bases 20x %"<80)], pos=2, srt=270, cex=.6, col="hotpink")
}

plot(samp$"Zero Coverage Targets %", main="Zero Coverage", ylab="% Targets with zero coverage", log="y", xlab="Sample", ylim=c(0.01,100))
abline(h=3, lty=2)
legend("bottomright", lty=2, legend="3% Targets Zero Coverage")
lowcoverage<-c(lowcoverage,samp$Sample[which(samp$"Zero Coverage">3)])
if(length(which(samp$"Zero Coverage Targets %">3))>0){
text(which(samp$"Zero Coverage Targets %">3), samp$"Zero Coverage Targets %"[which(samp$"Zero Coverage Targets %">3)], labels=samp$Sample[which(samp$"Zero Coverage Targets %">3)], pos=2, srt=270, cex=.6, col="hotpink")
}

print("Coverage stats done")
nofp<-lane$"External ID"[which(is.na(lane$"FP LOD"))]

if(length(which(is.na(lane$"FP LOD")))< nrow(lane)){

plot(lane$"FP Confident Calls"~as.factor(lane$"External ID"), xlab="sample", ylab="Multiplex level # FP calls", main="Fingerprint Calls/Sample Instance", xaxt="n")
medians<-tapply(lane$"FP Confident Calls",lane$"External ID", median, na.rm=TRUE)
points(as.factor(dimnames(medians)[[1]]),medians,col="red", lwd=2)
legend("topleft", legend="Median across sample instances", pch=1, lwd=2, col="red", lty=0)
poorFPcov<-dimnames(medians)[[1]][which(medians<5 )]
if(length(poorFPcov)>0){
text(which(medians<5), medians[which(medians<5)],poorFPcov, pos=2, srt=270, cex=.6, col="hotpink")
}

print("1 fp plot")
plot(100*(lane$"FP Confident Matching SNPs"/lane$"FP Confident Calls")~as.factor(lane$"External ID"), xlab="sample", ylab="Multiplex level % matching FP calls", main="% Confident calls matching for samples with low confident calls", xaxt="n", ylim=c(0,110))

print("2 fp plot")

plot(lane$"FP LOD"~as.factor(lane$"External ID"), xlab="sample", ylab="Sample Fingerprint LOD", main="Fingerprint Pass:Samples", xaxt="n")
offsamps<-lane$"External ID"[which(lane$"FP LOD"<(-3))]
lowfpLOD<-lane$"External ID"[which(lane$"FP LOD"<6)]

if(length(lowfpLOD)>0){
text(which(lane$"External ID" %in% lowfpLOD), lane$"FP_LOD"[which(lane$"FP LOD"<6)], labels=lowfpLOD, pos=2, srt=270, cex=.6, col="hotpink")
}
print("3 fp plot")

if(length(lowfpLOD)>0){
plot((lane$"FP Confident Calls"-lane$"FP Confident Matching SNPs")~as.factor(lane$"External ID"), main="Calls vs Matching Calls for Samples failing FP QC", ylab="# Mismatches", xlab="")
}
if(length(lowfpLOD)>0){
text(which(lane$"FP LOD"<6), lane$"FP_LOD"[which(lane$"FP LOD"<6)], labels=lowfpLOD, pos=2, srt=270, cex=.6, col="RED")
}


}else{
offsamps<-"NO FPDATA"
lowfpLOD<-"NO FP DATA"
poorFPcov<-"NO FP DATA"
}
print("FP stats done")

boxplot(samp$"Total SNPs", compsamp$"Total SNPs", names=c("Current Set", "All Sets"), ylab="Total SNPs per sample", main="Total SNPs")
standardQuants<-boxplot.stats(compsamp$"Total SNPs")$stats
offSNPs<-samp$Sample[which(samp$"Total SNPs" <standardQuants[1])]
offSNPs<-c(offSNPs, samp$Sample[which(samp$"Total SNPs" >standardQuants[5])])
if(length(offSNPs >0)){
  text(1, samp$"Total SNPs"[which(samp$Sample %in% offSNPs)], labels=offSNPs, pos=2, col="hot pink")
}
print("SNP stats done")

boxplot(samp$"dbSNP %", compsamp$"dbSNP %", names=c("Current Set", "All Sets"), ylab="% SNPs in dbSNP per sample", main="dbSNP Percentage")
standardQuants<-boxplot.stats(compsamp$"dbSNP %")$stats
offdbSNP<-samp$Sample[which(samp$"dbSNP %" <standardQuants[1])]
offdbSNP<-c(offdbSNP, samp$Sample[which(samp$"dbSNP %" >standardQuants[5])])
if(length(offdbSNP >0)){
  text(1, samp$"dbSNP %"[which(samp$Sample %in% offdbSNP)], labels=offdbSNP, pos=2, col="hot pink")
}
print("DBSNP stats done")

sampDuplication<-sub(pattern="Catch-.*: ", "",samp$"Library Duplication %")
sampDuplication<-as.numeric(sub("%", "", sampDuplication))
compsampDuplication<-sub(pattern="Catch-.*: ", "",compsamp$"Library Duplication %")
compsampDuplication<-as.numeric(sub("%", "", compsampDuplication))

boxplot(sampDuplication, compsampDuplication, names=c("Current Set", "All Sets"), ylab="% Duplication", main="Library Duplication")
standardQuants<-boxplot.stats(compsampDuplication)$stats
offDup<-samp$Sample[which(sampDuplication <standardQuants[1])]
offDup<-c(offDup, samp$Sample[which(sampDuplication >standardQuants[5])])
if(length(offDup >0)){
  text(1, sampDuplication[which(samp$Sample %in% offDup)], labels=offDup, pos=2, col="hot pink")
}
print("Duplication stats done")

allproblemsamples<-unique(c(lowcoverage, poorFPcov, offsamps, lowfpLOD, offSNPs, offdbSNP, offDup, missingLane, missingSamp))
problemMat<-matrix(c(rep("PASS", length(allproblemsamples)*9)), nrow=length(allproblemsamples))
rownames(problemMat)<-allproblemsamples
colnames(problemMat)<-c("low coverage", "low fp cov", "Identity Fail", "low FP LOD", "weird SNP count", "weird dbSNP %", "Duplicated", "Missing lane data", "missing agg data")
problemMat[which(rownames(problemMat) %in% lowcoverage),1]<-"FAIL"
problemMat[which(rownames(problemMat) %in% poorFPcov),2]<-"FAIL"
problemMat[which(rownames(problemMat) %in% offsamps),2]<-"FAIL"
problemMat[which(rownames(problemMat) %in% lowfpLOD),4]<-"FAIL"
problemMat[which(rownames(problemMat) %in% offSNPs),5]<-"FAIL"
problemMat[which(rownames(problemMat) %in% offdbSNP),6]<-"FAIL"
problemMat[which(rownames(problemMat) %in% offDup),7]<-"FAIL"
problemMat[which(rownames(problemMat) %in% missingLane),8]<-"FAIL"
problemMat[which(rownames(problemMat) %in% missingSamp),9]<-"FAIL"

textplot(problemMat, cex=.5)

write.table(problemMat, file=paste(cmdargs$qcout,"qc.table",sep="."), quote=FALSE, sep="\t")
print("no fp")
print(unique(nofp))



dev.off()
print("All stats done")
