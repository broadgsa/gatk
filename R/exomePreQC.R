data <- read.table('GoT2D_exomes_batch_005.tsv',header=T)

fingerprint_lods = list()
for(i in 1:nrow(data)) {
  fingerprint_lods[[as.character(data$sample[i])]] <- eval(parse(text=data$FINGERPRINT_LODS[i]))
}

fingerprint_lod_order = order(unlist(lapply(fingerprint_lods,median),use.names=F))

boxplot(fingerprint_lods[fingerprint_lod_order],las=3,main='Fingerprint LOD Scores By Sample',xlab='Sample',ylab='LOD Score Distribution',cex.axis=0.65)
