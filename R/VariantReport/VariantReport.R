source(paste(Sys.getenv("STING_DIR"), "/R/gsacommons.R", sep=""));

if (!interactive()) {
    options(warn=-1);

    library(getopt);

    spec = c(
        'evalRoot', 'e', 1, 'character', "root of the VariantEval files",
        'plotRoot', 'p', 1, 'character', "output root of the PDF file"
    );

    opt = getopt(matrix(spec, byrow=T, nrow=length(spec)/5));
} else {
    opt = list(
        evalRoot = "results/v9/with_1kg/recalibrated.with1KGSites.vcf.intermediate/merged.vcf.eval",
        plotRoot = "results/v9/with_1kg/recalibrated.with1KGSites.vcf.intermediate/merged.vcf.eval/plot"
    );
}

eval = read.eval(opt$evalRoot);

plot.begin(opt$plotRoot, "variantReport");

if (length(eval$CallsetNames) > 0) {
    # Venn diagram
    plot.callsetConcordance(eval);

    # Venn by AC
    plot.callsetConcordanceByAC(eval, novelty_name="all");
    plot.callsetConcordanceByAC(eval, novelty_name="known");
    plot.callsetConcordanceByAC(eval, novelty_name="novel");

    # Allele count spectrum
    plot.alleleCountSpectrum(eval, novelty_name="all");
    plot.alleleCountSpectrum(eval, novelty_name="known");
    plot.alleleCountSpectrum(eval, novelty_name="novel");

    # Ti/Tv spectrum
    plot.titvSpectrum(eval, novelty_name="all");
    plot.titvSpectrum(eval, novelty_name="known");
    plot.titvSpectrum(eval, novelty_name="novel");

    # Per-sample
    #plot.variantsPerSample(eval);
} else {
    plot.variantsPerSample(eval, novelty_name="all");
    plot.variantsPerSample(eval, novelty_name="known");
    plot.variantsPerSample(eval, novelty_name="novel");
}

plot.variantsPerSample(eval, novelty_name="all");
plot.variantsPerSample(eval, novelty_name="known");
plot.variantsPerSample(eval, novelty_name="novel");

allVariants = eval$CountVariants[which(eval$CountVariants$jexl_expression == "none" & eval$CountVariants$filter_name == "called" & eval$CountVariants$novelty_name == "all"),]$nVariantLoci;
knownVariants = eval$CountVariants[which(eval$CountVariants$jexl_expression == "none" & eval$CountVariants$filter_name == "called" & eval$CountVariants$novelty_name == "known"),]$nVariantLoci;
novelVariants = eval$CountVariants[which(eval$CountVariants$jexl_expression == "none" & eval$CountVariants$filter_name == "called" & eval$CountVariants$novelty_name == "novel"),]$nVariantLoci;

allTiTv = eval$TiTv[which(eval$TiTv$jexl_expression == "none" & eval$TiTv$filter_name == "called" & eval$TiTv$novelty_name == "all"),]$ti.tv_ratio;
knownTiTv = eval$TiTv[which(eval$TiTv$jexl_expression == "none" & eval$TiTv$filter_name == "called" & eval$TiTv$novelty_name == "known"),]$ti.tv_ratio;
novelTiTv = eval$TiTv[which(eval$TiTv$jexl_expression == "none" & eval$TiTv$filter_name == "called" & eval$TiTv$novelty_name == "novel"),]$ti.tv_ratio;

suppressPackageStartupMessages(library(gplots));

variantsummary = as.data.frame(t(matrix(c(allVariants, allTiTv, knownVariants, knownTiTv, novelVariants, novelTiTv), nrow=2)));
colnames(variantsummary) = c("counts", "ti/tv");
rownames(variantsummary) = c("all", "known", "novel");

textplot(variantsummary);

plot.end(opt$plotRoot);
