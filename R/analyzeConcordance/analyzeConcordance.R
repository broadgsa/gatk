#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

base_name = args[1]
input = args[2]

d <- read.table(input, header=T)
# separate the data into filtered and unfiltered

d.filtered <- d[d$filter_type=="filtered",]
d.unfiltered <- d[d$filter_type=="unfiltered",]

if (nrow(d.filtered) > 0) {
	d.display <- d.filtered
} else {
	d.display <- d.unfiltered
}

#
# Plot histograms of the known versus novel Ti/Tv
#

outfile = paste(base_name, ".histograms.png", sep="")

if (nrow(d.filtered) > 0) {
  nFilterTypes <- 2
} else {
  nFilterTypes <- 1
}

bitmap(outfile, width=600, height=(300 * nFilterTypes), units="px")
par(cex=1.1, mfrow=c(1 * nFilterTypes,2))
nbreaks <- 20
color <- "grey"
xlim <- c(0,4)

hist(d.unfiltered$known_titv, nbreaks, col=color, xlim=xlim)
hist(d.unfiltered$novel_titv, nbreaks, col=color, xlim=xlim)

if (nrow(d.filtered) > 0) {
  hist(d.filtered$known_titv, nbreaks, col=color, xlim=xlim)
  hist(d.filtered$novel_titv, nbreaks, col=color, xlim=xlim)
}

dev.off()

#
# Plot samples in order of novel Ti/Tv versus known Ti/Tv
#

outfile = paste(base_name, ".novel_vs_known_titv.png", sep="")

bitmap(outfile, width=600, height=600, units="px")

d.display <- d.display[order(d.display$novel_titv),]
plot(1:length(d.display$known_titv),d.display$known_titv,type="b",col="blue",ylim=c(0,4), xlab="Sample #", ylab="Ti / Tv")
points(1:length(d.display$novel_titv),d.display$novel_titv,type="b",col="red",ylim=c(0,4))
legend("bottomright", c("known","novel"), col=c("blue","red"), pch=21)

dev.off()
