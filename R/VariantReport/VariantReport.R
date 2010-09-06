source(paste(Sys.getenv("STING_DIR"), "/R/gsacommons.R", sep=""));

if (interactive()) {
    if (!exists("plotRoot")) {
        plotRoot = "test.plot";
    }
} else {
    args = commandArgs(TRUE);

    evalRoot = args[1];
    plotRoot = args[2];
}

eval = read.eval(evalRoot);

plot.begin(plotRoot, "variantReport");

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

plot.end(plotRoot);
