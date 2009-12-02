#!/broad/tools/apps/R-2.6.0/bin/Rscript

args <- commandArgs(TRUE)
verbose = TRUE

input = args[1]

#outfile = paste(input, ".qual_diff_v_dinuc.png", sep="")
#png(outfile, height=7, width=7, units="in", res=72) #height=1000, width=680)
outfile = paste(input, ".qual_diff_v_dinuc.pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
#in_dinuc = paste(input, ".quality_difference_v_dinucleotide.csv", sep="")
#d <- read.csv(input)
d <- read.table(input, header=T)
plot(d$Dinuc, d$Qempirical_Qreported, type="l", ylab="Empirical - Reported Quality", xlab="Dinucleotide", ylim=c(-10,10))
