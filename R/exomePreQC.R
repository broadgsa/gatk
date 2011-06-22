args = commandArgs(TRUE)
onCMDLine = ! is.na(args[1])

if ( onCMDLine ) {
  reference_dataset = '/Users/mhanna/metrics.perSample.formatted.table'
  inputTSV = args[1]
  outputPDF = args[2]
} else {
  reference_dataset = '/Users/mhanna/metrics.perSample.formatted.table'
  inputTSV = 'GoT2D_exomes_batch_005.tsv'
  outputPDF = 'T2D.pdf'
}

require('ggplot2')

data <- read.table(inputTSV,header=T)

complete <- read.table(reference_dataset,header=T)
novel <- subset(complete,exon_intervals == "whole_exome_agilent_1.1_refseq_plus_3_boosters"&Novelty=="novel"&FunctionalClass=="all")
selected_samples <- novel$Sample %in% data$sample
novel_with_highlights <- cbind(novel,selected_samples)

if(onCMDLine) {
    fingerprint_lods = list()
    for(i in 1:nrow(data)) {
      fingerprint_lods[[as.character(data$sample[i])]] <- eval(parse(text=data$FINGERPRINT_LODS[i]))
    }

    fingerprint_lod_order = order(unlist(lapply(fingerprint_lods,median),use.names=F))

    pdf(outputPDF)
    boxplot(fingerprint_lods[fingerprint_lod_order],las=3,main='Fingerprint LOD Scores By Sample',xlab='Sample',ylab='LOD Score Distribution',cex.axis=0.65)

    qplot(Sample,Selected_Bases_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='On+Near Bait Bases/PF Bases Aligned per Sample')
    qplot(Sample,Mean_Target_Coverage,data=novel_with_highlights,color=selected_samples) + opts(title='Mean Target Coverage per Sample')
    qplot(Sample,Zero_Coverage_Targets_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% of Targets with <2x Coverage per Sample')
    qplot(Sample,Fold_80_Base_Penalty,data=novel_with_highlights,color=selected_samples) + opts(title='Fold 80 Base Penalty per Sample')
    qplot(Sample,Target_Bases_20x_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% Target Bases Achieving >20x Coverage per Sample')
    qplot(Sample,PF_Reads_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% PF Reads Aligned per Sample')
    qplot(Sample,PF_HQ_Error_Rate,data=novel_with_highlights,color=selected_samples) + opts(title='% HQ Bases mismatching the Reference per Sample')
    qplot(Sample,Mean_Read_Length,data=novel_with_highlights,color=selected_samples) + opts(title='Median Read Length per Sample')
    qplot(Sample,Bad_Cycles,data=novel_with_highlights,color=selected_samples) + opts(title='# Bad Cycles per Sample')
    qplot(Sample,Strand_Balance_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% PF Reads Aligned to the + Strand per Sample')
    qplot(Sample,Total_SNPs,data=novel_with_highlights,color=selected_samples) + opts(title='# SNPs called per Sample')
    qplot(Sample,dbSNP_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% SNPs in dbSNP per Sample')
    qplot(PCT_DBSNP,data=data,geom="histogram") + opts(title='% SNPs in dbSNP per Sample')
    dev.off()
} else {
     print('Plotting command-line arguments')
     qplot(Sample,PF_Reads_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% PF Reads Aligned per Sample')
}

#qplot(Sample,Library_Size_HS,data=novel_with_highlights,color=selected_samples) + opts(title='Hybrid Sequencing Library Size per Sample')
#qplot(Sample,MEDIAN_INSERT_SIZE,data=novel_with_highlights,color=selected_samples) + opts(title='Median Insert Size per Sample')
#qplot(Sample,PCT_CHIMERAS,data=novel_with_highlights,color=selected_samples) + opts(title='% Chimera Read Pairs per Sample')
#qplot(Sample,PCT_ADAPTER,data=novel_with_highlights,color=selected_samples) + opts(title='% Unaligned Reads Matching an Adapter Sequence per Sample')
#qplot(Sample,NOVEL_SNPS,data=novel_with_highlights,color=selected_samples) + opts(title='# Novel SNPs called per Sample')
#qplot(Sample,DBSNP_TITV,data=novel_with_highlights,color=selected_samples) + opts(title='TiTv of SNPs in dbSNP per Sample')
