args = commandArgs(TRUE)
onCMDLine = ! is.na(args[1])

if ( onCMDLine ) {
  inputTSV = args[1]
  outputPDF = args[2]
}

require('ggplot2')

inputTSV = "GoT2D_exomes_batch_005.tsv"
data <- read.table(inputTSV,header=T)

fingerprint_lods = list()
for(i in 1:nrow(data)) {
  fingerprint_lods[[as.character(data$sample[i])]] <- eval(parse(text=data$FINGERPRINT_LODS[i]))
}

fingerprint_lod_order = order(unlist(lapply(fingerprint_lods,median),use.names=F))

pdf(outputPDF)
boxplot(fingerprint_lods[fingerprint_lod_order],las=3,main='Fingerprint LOD Scores By Sample',xlab='Sample',ylab='LOD Score Distribution',cex.axis=0.65)

complete <- read.table('/Users/mhanna/metrics.perSample.formatted.table',header=T)
novel <- subset(complete,exon_intervals == "whole_exome_agilent_1.1_refseq_plus_3_boosters"&Novelty=="novel"&FunctionalClass=="all")
selected_samples <- novel$Sample %in% data$sample
novel_with_highlights <- cbind(novel,selected_samples)

qplot(Sample,Selected_Bases_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='On+Near Bait Bases/PF Bases Aligned per Sample')
qplot(PCT_SELECTED_BASES,data=data,geom="histogram") + opts(title='On+Near Bait Bases (Distribution)')
qplot(Sample,Mean_Target_Coverage,data=novel_with_highlights,color=selected_samples) + opts(title='Mean Target Coverage per Sample')
qplot(MEAN_TARGET_COVERAGE,data=data,geom="histogram") + opts(title='Mean Target Coverage (Distribution)')
qplot(Sample,Zero_Coverage_Targets_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% of Targets with <2x Coverage per Sample')
qplot(ZERO_CVG_TARGETS_PCT,data=data,geom="histogram") + opts(title='% of Targets with <2x Coverage (Distribution)')
qplot(Sample,Fold_80_Base_Penalty,data=novel_with_highlights,color=selected_samples) + opts(title='Fold 80 Base Penalty per Sample')
qplot(FOLD_80_BASE_PENALTY,data=data,geom="histogram") + opts(title='Fold 80 Base Penalty (Distribution)')
qplot(Sample,Target_Bases_2x_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% Target Bases Achieving >2x Coverage per Sample')
qplot(PCT_TARGET_BASES_2X,data=data,geom="histogram") + opts(title='% Target Bases Achieving >2x Coverage (Distribution)')
qplot(Sample,Target_Bases_10x_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% Target Bases Achieving >10x Coverage per Sample')
qplot(PCT_TARGET_BASES_10X,data=data,geom="histogram") + opts(title='% Target Bases Achieving >10x Coverage (Distribution)')
qplot(Sample,Target_Bases_20x_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% Target Bases Achieving >20x Coverage per Sample')
qplot(PCT_TARGET_BASES_20X,data=data,geom="histogram") + opts(title='% Target Bases Achieving >20x Coverage (Distribution)')
qplot(Sample,Target_Bases_30x_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% Target Bases Achieving >30x Coverage per Sample')
qplot(PCT_TARGET_BASES_30X,data=data,geom="histogram") + opts(title='% Target Bases Achieving >30x Coverage (Distribution)')
qplot(Sample,PF_Reads_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% PF Reads Aligned per Sample')
qplot(PCT_PF_READS_ALIGNED,data=data,geom="histogram") + opts(title='% PF Reads Aligned (Distribution)')
qplot(Sample,PF_HQ_Error_Rate,data=novel_with_highlights,color=selected_samples) + opts(title='% HQ Bases mismatching the Reference per Sample')
qplot(PF_HQ_ERROR_RATE,data=data,geom="histogram") + opts(title='% HQ Bases mismatching the Reference (Distribution)')
qplot(Sample,Mean_Read_Length,data=novel_with_highlights,color=selected_samples) + opts(title='Median Read Length per Sample')
qplot(MEAN_READ_LENGTH,data=data,geom="histogram") + opts(title='Median Read Length (Distribution')
qplot(Sample,Bad_Cycles,data=novel_with_highlights,color=selected_samples) + opts(title='# Bad Cycles per Sample')
qplot(BAD_CYCLES,data=data,geom="histogram") + opts(title='# Bad Cycles (Distribution)')
qplot(Sample,Strand_Balance_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% PF Reads Aligned to the + Strand per Sample')
qplot(STRAND_BALANCE,data=data,geom="histogram") + opts(title='% PF Reads Aligned to the + Strand (Distribution)')
qplot(Sample,Total_SNPs,data=novel_with_highlights,color=selected_samples) + opts(title='# SNPs called per Sample')
qplot(TOTAL_SNPS,data=data,geom="histogram") + opts(title='# SNPs called (Distribution)')
qplot(Sample,dbSNP_Pct,data=novel_with_highlights,color=selected_samples) + opts(title='% SNPs in dbSNP per Sample')
qplot(PCT_DBSNP,data=data,geom="histogram") + opts(title='% SNPs in dbSNP per Sample')
dev.off()

#qplot(Sample,Library_Size_HS,data=novel_with_highlights,color=selected_samples) + opts(title='Hybrid Sequencing Library Size per Sample')
qplot(HS_LIBRARY_SIZE,data=data) + opts(title='Hybrid Sequencing Library Size (Distribution)')
#qplot(Sample,MEDIAN_INSERT_SIZE,data=novel_with_highlights,color=selected_samples) + opts(title='Median Insert Size per Sample')
qplot(MEDIAN_INSERT_SIZE,data=data) + opts(title='Median Insert Size (Distribution)')
#qplot(Sample,PCT_CHIMERAS,data=novel_with_highlights,color=selected_samples) + opts(title='% Chimera Read Pairs per Sample')
qplot(PCT_CHIMERAS,data=data) + opts(title='% Chimera Read Pairs (Distribution)')
#qplot(Sample,PCT_ADAPTER,data=novel_with_highlights,color=selected_samples) + opts(title='% Unaligned Reads Matching an Adapter Sequence per Sample')
qplot(PCT_ADAPTER,data=data) + opts(title='% Unaligned Reads Matching an Adapter Sequence (Distribution)')
#qplot(Sample,NOVEL_SNPS,data=novel_with_highlights,color=selected_samples) + opts(title='# Novel SNPs called per Sample')
qplot(NOVEL_SNPS,data=data) + opts(title='# Novel SNPs called (Distribution)')
#qplot(Sample,DBSNP_TITV,data=novel_with_highlights,color=selected_samples) + opts(title='TiTv of SNPs in dbSNP per Sample')
qplot(DBSNP_TITV,data=data) + opts(title='TiTv of SNPs in dbSNP (Distribution)')
