source(paste(Sys.getenv("STING_DIR"), "/R/gsacommons.R", sep=""));

args = commandArgs(TRUE);

evalRoot = args[1];
#if (is.na(evalRoot)) { evalRoot = "/home/radon01/kiran/scr1/projects/ESPGabrielDownsampling/results/v1/merged.vcf.eval120/eval"; }
if (is.na(evalRoot)) { evalRoot = "/home/radon01/kiran/scr1/projects/ESPGabrielDownsampling/results/v1/merged.vcf.eval.v2/eval"; }

plotRoot = args[2];
if (is.na(plotRoot)) { plotRoot = "plot"; }

eval = read.eval(evalRoot);

# Venn diagram
plot.begin(plotRoot, "venn");
plot.callsetConcordance(eval);
plot.end(plotRoot);

# Venn by AC
plot.begin(plotRoot, "venn_by_ac.all", width=12, height=8);
plot.callsetConcordanceByAC(eval, novelty_name="all");
plot.end(plotRoot);

plot.begin(plotRoot, "venn_by_ac.known", width=12, height=8);
plot.callsetConcordanceByAC(eval, novelty_name="known");
plot.end(plotRoot);

plot.begin(plotRoot, "venn_by_ac.novel", width=12, height=8);
plot.callsetConcordanceByAC(eval, novelty_name="novel");
plot.end(plotRoot);

# Allele count spectrum
plot.begin(plotRoot, "acs.all", width=12, height=8);
plot.alleleCountSpectrum(eval, novelty_name="all");
plot.end(plotRoot);

plot.begin(plotRoot, "acs.known", width=12, height=8);
plot.alleleCountSpectrum(eval, novelty_name="known");
plot.end(plotRoot);

plot.begin(plotRoot, "acs.novel", width=12, height=8);
plot.alleleCountSpectrum(eval, novelty_name="novel");
plot.end(plotRoot);

# Ti/Tv spectrum
plot.begin(plotRoot, "titv.all", width=12, height=8);
plot.titvSpectrum(eval, novelty_name="all");
plot.end(plotRoot);

plot.begin(plotRoot, "titv.known", width=12, height=8);
plot.titvSpectrum(eval, novelty_name="known");
plot.end(plotRoot);

plot.begin(plotRoot, "titv.novel", width=12, height=8);
plot.titvSpectrum(eval, novelty_name="novel");
plot.end(plotRoot);

# Per-sample
plot.begin(plotRoot, "variants_per_sample", width=12, height=8);
plot.variantsPerSample(eval);
plot.end(plotRoot);
