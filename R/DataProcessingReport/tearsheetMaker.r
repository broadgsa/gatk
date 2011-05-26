source("/humgen/gsa-pipeline/.repository/R/DataProcessingReport/Tearsheet.R")
cmdargs = gsa.getargs(
    list(
        title = list(value=NA, doc="Title for the tearsheet"),
        tsv = list(value=NA, doc="pipeline tsv file"),
        evalroot = list(value=NA, doc="VariantEval file base (everything before the .eval)"),
        tearout = list(value=NA, doc="Output path for tearsheet PDF")#,
    ),
    doc="Creates a tearsheet"
);

read.delim(cmdargs$tsv, header=FALSE)->settable

squids<-unique(settable[,1])

lane<-data.frame()
samp<-data.frame()
for(squid in squids){
  gsa.read.squidmetrics(squid, TRUE)->lanemetrics 
  addlanes<-lanemetrics[which(lanemetrics$"External ID" %in% settable[,2]),]
  gsa.read.squidmetrics(squid, FALSE)->samplemetrics 
  addsamps<-samplemetrics[which(samplemetrics$"Sample" %in% settable[,2]),]
  lane<-rbind(lane, addlanes)
  samp<-rbind(samp, addsamps)
}
print("Picard Data Obtained...")
gsa.read.gatkreport(paste(cmdargs$evalroot, ".eval", sep=""))->basiceval 
gsa.read.gatkreport(paste(cmdargs$evalroot, ".extraFC.eval", sep=""))->FCeval 
gsa.read.gatkreport(paste(cmdargs$evalroot, ".extraSA.eval", sep=""))->SAeval 
print("Evals read")
  pdf(file= cmdargs$tearout, width=22, height=17, pagecentre=TRUE, pointsize=24) 
    print("PDF created...")
tearsheet()
dev.off()