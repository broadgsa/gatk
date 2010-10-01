source(paste(Sys.getenv("STING_DIR"), "/R/gsacommons.R", sep=""));

if (!interactive()) {
    options(warn=-1);

    library(getopt);

    spec = c(
        'evalRoot', 'e', 1, 'character', "root of the VariantEval files",
        'plotRoot', 'p', 1, 'character', "output root of the PDF file",
        'title',    't', 1, 'character', "title of the report",
        'author',   'a', 1, 'character', "author that should be listed on the report"
    );

    opt = getopt(matrix(spec, byrow=T, nrow=length(spec)/5));
} else {
    opt = list(
        evalRoot = "results/v9/with_1kg/recalibrated.with1KGSites.vcf.intermediate/merged.vcf.eval",
        plotRoot = "results/v9/with_1kg/recalibrated.with1KGSites.vcf.intermediate/merged.vcf.eval/plot",
        title = "Test"
    );
}

eval = read.eval(opt$evalRoot);

plot.begin(opt$plotRoot, "variantReport");

plot.titlePage(opt$title, opt$author);

plot.variantTable(eval);

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
    #plot.variantsPerSample(eval, novelty_name="all");
    #plot.variantsPerSample(eval, novelty_name="known");
    #plot.variantsPerSample(eval, novelty_name="novel");
}

plot.end(opt$plotRoot);
