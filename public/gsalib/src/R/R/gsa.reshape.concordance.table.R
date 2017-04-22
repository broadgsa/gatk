gsa.reshape.concordance.table <- function(report, table.name="GenotypeConcordance_Counts", sample.name="ALL") {
  if (!is.null(table.name)) {
    data <- report[[table.name]]
  }
  if (is.null(table.name)) {
    data <- report
  }
  d <- data[data$Sample==sample.name,2:(length(data[1,])-1)]
  
  possible.genotypes <- c('NO_CALL', 'HOM_REF', 'HET', 'HOM_VAR', 'UNAVAILABLE', 'MIXED')
  combinations <- outer(possible.genotypes, possible.genotypes, function(a,b) {paste(a,b,sep='_')})
  existing.combi <- matrix(combinations %in% colnames(d), nrow=length(possible.genotypes))
  eval.genotypes <- apply(existing.combi, 1, any)
  comp.genotypes <- apply(existing.combi, 2, any)
  
  m <- matrix(d, nrow=sum(eval.genotypes), byrow=T)
  dimnames(m) <- list(EvalGenotypes=possible.genotypes[eval.genotypes],
                      CompGenotypes=possible.genotypes[comp.genotypes])
  m
}
