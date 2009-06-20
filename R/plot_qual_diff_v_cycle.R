#!/broad/tools/apps/R-2.6.0/bin/Rscript

args <- commandArgs(TRUE)
verbose = TRUE

input = args[1]

#X11(width=7, height=14)
#outfile = paste(input, ".qual_diff_v_cycle.png", sep="")
#png(outfile, height=7, width=7, units="in", res=72) #height=1000, width=680)
outfile = paste(input, ".qual_diff_v_cycle.pdf", sep="")
pdf(outfile, height=7, width=7)
par(cex=1.1)
c <- read.table(input, header=T)
plot(c$Cycle, c$Qempirical_Qreported, type="l", ylab="Empirical - Reported Quality", xlab="Cycle", ylim=c(-10, 10))

