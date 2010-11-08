#d <- read.table("sim_calls.table", header=T)

d$sim.VAR <- d$sim.AC > 0
d$called.VAR <- d$called.AC > 0

QS = unique(d$sim.Q)
MODES = unique(d$sim.MODE)
NS = unique(d$called.AN / 2)
DEPTHS = unique(d$sim.DP)

results <- expand.grid(Q = QS, mode = MODES, nSamples = NS, depth = DEPTHS)
results$sensitivity = 0
results$specificity = 0

determineRates <- function(raw, Q, mode, depth) {
    sub <- subset(raw, sim.Q == Q & sim.MODE == mode & sim.DP == depth)
    ct <- table(sub$called.VAR, sub$sim.VAR, dnn = c("called.VAR", "sim.VAR"))
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
    boxplot(called.AC ~ sim.AC, data = subset(d, called.DP == depth * NS), main = paste("Depth of coverage ", depth), xlab = "Simulation AC", ylab = "Called AC")

print(results)

par(mfcol=c(2,1))
for ( Qt in QS ) {
    x <- subset(results, Q == Qt)
    print(x)
    plot(x$depth, x$sensitivity, type="b")
    plot(x$depth, x$specificity, type="b")
}
