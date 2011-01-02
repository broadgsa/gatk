args = commandArgs(TRUE);

RUNME = T
onCMDLine = ! is.na(args[1])
DATA_FILE = args[1]
DESCRIPTION = args[2]
OUTPUT_PDF = paste(DATA_FILE, ".pdf", sep="")

if ( onCMDLine ) { 
   print(paste("Reading data from", DATA_FILE))
   d = read.table(DATA_FILE, header=T)
} 

if ( onCMDLine ) pdf(OUTPUT_PDF)

generateOneReport <- function(d) {
    qs = quantile(d$processing.speed, probs = c(0.05, 0.5, 0.95))
    plot(d$elapsed.time, d$processing.speed, main=DESCRIPTION, xlab="Elapsed time (sec)", ylab="Processing speed (seconds per 1M units)", ylim=c(qs[1], qs[3]), type="b", col="cornflowerblue", lwd=2)
    abline(h=qs[2], lty=2)
}

if ( RUNME ) {
    generateOneReport(d)
}

if ( onCMDLine ) dev.off()



