source("/humgen/gsa-pipeline/.repository/R/DataProcessingReport/qcplots.r")
suppressMessages(library(gplots));
def.par <- par(no.readonly = TRUE)


cmdargs = gsa.getargs(
    list(
        tsv = list(value=NA, doc="pipeline tsv file"),
        evalroot = list(value=NA, doc="VariantEval file base (everything before the .eval)"),
        reportout = list(value=NA, doc="Output path for report PDF")#,
    ),
    doc="Creates a variant report"
);

read.delim(cmdargs$tsv, header=FALSE)->settable

squids<-unique(settable[,1])


gsa.read.gatkreport(paste(cmdargs$evalroot, ".eval", sep=""))->basiceval 
gsa.read.gatkreport(paste(cmdargs$evalroot, ".extraSA.eval", sep=""))->SAeval 
print("Evals read")

pdf(file= cmdargs$reportout, width=22, height=17, pagecentre=TRUE, pointsize=24) 
    print("PDF created...")


path="."
weirdos<-which(SAeval$TiTvVariantEvaluator$Sample %in% SAeval$TiTvVariantEvaluator$Sample[which(SAeval$TiTvVariantEvaluator$tiTvRatio <2)])

novelAC(SAeval)
knownAC(SAeval)
AllAC(SAeval)
layout(matrix(c(6,1, 2,3, 4, 5), nrow=6), heights=c(1, 1, 1, 1, 1,1)) 
textplot("Sample Novel TiTv ranges should be above 2, as they are in previous datasets. \nSamples with lower TiTv data are flagged in subsequent plots with hot pink labels, and listed below:")
textplot(paste(unique(SAeval$TiTvVariantEvaluator$Sample[weirdos]), collapse=", "), halign="left")
textplot("Problem Samples frequently have unusually high or low numbers of variants.")
textplot("Samples with unusually high numbers of novel variants may be from different populations, and, as such, should have higher heterozygosity. \nIf this is not the case, there may be problems with the samples.")
textplot("Unusually high numbers of variants with low allele counts may indicate variants generated from problematic samples.")
textplot("Notes for interpreting QC data:")
dev.off()
