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
d.good <- c[c$nMismatches >= 100,]
d.100 <- c[c$nMismatches < 100,]
plot(d.good$Cycle, d.good$Qempirical_Qreported, type="l", ylab="Empirical - Reported Quality", xlab="Cycle", col="blue", ylim=c(-10, 10))
points(d.100$Cycle, d.100$Qempirical_Qreported, type="p", col="lightblue", pch=3)
#points(d.1000$Cycle, d.1000$Qempirical_Qreported, type="p", col="cornflowerblue", pch=16)

