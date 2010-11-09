require("plotrix")
args = commandArgs(TRUE);

onCMDLine = ! is.na(args[1])

file = "sim_calls.table"
info = "interactive R"
if ( onCMDLine ) {
    file = args[1]
    d <- read.table(file, header=T)
    pdf(args[2])
    info = args[3]
}

d$sim.VAR <- d$sim.AC > 0
d$called.VAR <- d$called.AC > 0

QS = unique(d$sim.Q)
MODES = unique(d$sim.MODE)
NS = unique(d$called.AN / 2)
DEPTHS = unique(d$sim.DP)

addSection <- function(name) {
    par("mar", c(5, 4, 4, 2))
    frame()
    title(name, cex=2)
}   

addSection(paste("Calling performance report: nSamples = ", NS, "\n  info:", info))

results <- expand.grid(Q = QS, mode = MODES, nSamples = NS, depth = DEPTHS)
results$sensitivity = 0
results$specificity = 0

determineRates <- function(raw, Q, mode, depth) {
    sub <- subset(raw, sim.Q == Q & sim.MODE == mode & sim.DP == depth)
    print(c(Q,mode,depth, dim(sub)))
    ct <- table(sub$called.VAR, sub$sim.VAR, dnn = c("called.VAR", "sim.VAR"), useNA = "always")
    print(ct)
    sensitivity = ct[2,2] / sum(ct[,2]) 
    specificity = ct[1,1] / sum(ct[,1])
    list(sensitivity = sensitivity, specificity = specificity, ct = ct)
}

for ( i in 1:(dim(results)[1]) ) {
    r <- results[i,]
    x <- determineRates(d, r$Q, r$mode, r$depth)
    results[i,]$sensitivity = x$sensitivity
    results[i,]$specificity = x$specificity
}

for ( depth in DEPTHS ) 
    boxplot(called.AC ~ sim.AC, data = subset(d, called.DP == depth * NS), main = paste("Depth of coverage ", depth), xlab = "Simulation AC", ylab = "Called AC", outwex=0.5, col = "cornflowerblue")

print(results)

par(mfcol=c(2,1))
for ( Qt in QS ) {
    x <- subset(results, Q == Qt)
    print(x)
    plot(x$depth, x$sensitivity, type="b", main = paste("Q score", Qt), xlab = "Depth", ylab="Sensitivity")
    plot(x$depth, x$specificity, type="b", xlab = "Depth", ylab="Specificity")
}

par(mfcol=c(1,1))
plot(0,0, type="n", frame.plot=F, ann=F, axes=F)
addtable2plot(-1, -1, data.frame(Q=results$Q, mode=results$mode, depth=results$depth, sensitivity=format(results$sensitivity, digits=2), specificity = format(results$specificity, digits=2)))


if ( onCMDLine ) dev.off()

