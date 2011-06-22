require('ggplot2')

data <- read.table('GoT2D_exomes_batch_005.tsv',header=T)

fingerprint_lods = list()
for(i in 1:nrow(data)) {
  fingerprint_lods[[as.character(data$sample[i])]] <- eval(parse(text=data$FINGERPRINT_LODS[i]))
}

fingerprint_lod_order = order(unlist(lapply(fingerprint_lods,median),use.names=F))

pdf('T2D.pdf')
boxplot(fingerprint_lods[fingerprint_lod_order],las=3,main='Fingerprint LOD Scores By Sample',xlab='Sample',ylab='LOD Score Distribution',cex.axis=0.65)
qplot(sample,GENOME_SIZE,data=data) + opts(title='Genome Size per Sample')
qplot(sample,PCT_SELECTED_BASES,data=data) + opts(title='On+Near Bait Bases/PF Bases Aligned per Sample')
qplot(sample,MEAN_TARGET_COVERAGE,data=data) + opts(title='Mean Target Coverage per Sample')
qplot(sample,ZERO_CVG_TARGETS_PCT,data=data) + opts(title='% of Targets with <2x Coverage per Sample')
qplot(sample,FOLD_80_BASE_PENALTY,data=data) + opts(title='Fold 80 Base Penalty per Sample')
qplot(sample,HS_LIBRARY_SIZE,data=data) + opts(title='Hybrid Sequencing Library Size per Sample')
qplot(sample,PCT_PF_READS_ALIGNED,data=data) + opts(title='% PF Reads Aligned per Sample')
qplot(sample,PF_HQ_ERROR_RATE,data=data) + opts(title='% HQ Bases mismatching the Reference  per Sample')
qplot(sample,MEAN_READ_LENGTH,data=data) + opts(title='Median Read Length per Sample')
qplot(sample,MEDIAN_INSERT_SIZE,data=data) + opts(title='Median Insert Size per Sample')
qplot(sample,BAD_CYCLES,data=data) + opts(title='# Bad Cycles per Sample')
qplot(sample,STRAND_BALANCE,data=data) + opts(title='% PF Reads Aligned to the + Strand per Sample')
qplot(sample,PCT_CHIMERAS,data=data) + opts(title='% Chimera Read Pairs per Sample')
qplot(sample,PCT_ADAPTER,data=data) + opts(title='% Unaligned Reads Matching an Adapter Sequence per Sample')
qplot(sample,TOTAL_SNPS,data=data) + opts(title='# SNPs called per Sample')
qplot(sample,NOVEL_SNPS,data=data) + opts(title='# Novel SNPs called per Sample')
qplot(sample,PCT_DBSNP,data=data) + opts(title='% SNPs in dbSNP per Sample')
qplot(sample,DBSNP_TITV,data=data) + opts(title='TiTv of SNPs in dbSNP per Sample')
dev.off()

