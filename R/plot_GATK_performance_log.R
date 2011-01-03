args = commandArgs(TRUE);

RUNME = T
onCMDLine = ! is.na(args[1])
DATA_FILE = args[1]
DESCRIPTION = args[2]
#OUTPUT_PDF = paste(DATA_FILE, ".pdf", sep="")

MAX_POINTS = 100000

if ( onCMDLine ) { 
   print(paste("Reading data from", DATA_FILE))
   d = read.table(DATA_FILE, header=T)
} 

#if ( onCMDLine ) pdf(OUTPUT_PDF)

vec.margin <- function(x) {
    l = length(x)
    d = x[-1] - x[1:(l-1)]
    c(x[1], d[1:(l-1)])
}

everyNth <- function(x, n) {
    l = dim(x)[1]
    m = ceiling(l / n)
    print(m)
    keep = 1:l %% m == 0
    x[keep,]
}


l = length(d$units.processed)
d$units.processed.margin = vec.margin(d$units.processed)
#prev = 0
#for ( i in 1:l ) {
#    cur = d$units.processed[i]
#    d[i,]$units.processed.margin = cur - prev
#    prev = cur
#}

generateOneReport <- function(d) {
    qs = quantile(d$processing.speed, probs = c(0.01, 0.5, 0.99))

    # unit processing time
    if ( onCMDLine ) png(paste(DATA_FILE, ".speed.png", sep=""), width=1080, height=1080)
    dpoints = everyNth(d, MAX_POINTS)
    plot(dpoints$elapsed.time, dpoints$processing.speed, main=DESCRIPTION, xlab="Elapsed time (sec)", ylab="Processing speed (seconds per 1M units)", ylim=c(qs[1], qs[3]), type="b", col="cornflowerblue", lwd=2)
    abline(h=qs[2], lty=2)
    if ( onCMDLine ) dev.off()

    # instantaneous processing speed
    if ( onCMDLine ) png(paste(DATA_FILE, ".marginal.png", sep=""), width=1080, height=1080)
    running_median_window = 101
    rm = runmed(d$units.processed.margin, running_median_window)
    POINT_COL = "#0000AA33"
    plot(dpoints$elapsed.time, dpoints$units.processed.margin, main=DESCRIPTION, xlab="Elapsed time (sec)", ylab="Units processed in last timing interval", type="p", cex = 0.5, col=POINT_COL)
    lines(d$elapsed.time, rm, lwd=3, col="red")
    legend("topleft", c("Observations", "101-elt running median"), fill=c(POINT_COL, "red"))
    if ( onCMDLine ) dev.off()
}

if ( RUNME ) {
    generateOneReport(d)
}




