theta = 10^-3

Fm = 392
Fs = 44

L = 101

meanDepth = 65
nReadsToPhase = 1

params = paste("meanDepth= ", meanDepth, ", nReadsToPhase= ", nReadsToPhase, ", L= ", L, ", theta= ", theta, ", Fm= ", Fm, ", Fs= ", Fs, sep="")

#
#options(warn=2)
#options(error=recover)
#

MAX_WINDOW_SIZE = 10

WINDOW_SIZES = 2:MAX_WINDOW_SIZE
NUM_WINDOW_SIZES = length(WINDOW_SIZES)

pPhaseWindow = as.vector(matrix(data = -1, nrow = 1, ncol = NUM_WINDOW_SIZES))
for (i in 1:NUM_WINDOW_SIZES) {
	n = WINDOW_SIZES[i]
	phaseIndex = middleOfWindowIndex(n)
	pPhaseWindow[i] = pDirectlyPhaseHetPairUsingWindow(meanDepth, nReadsToPhase, L, theta, Fm, Fs, n, phaseIndex)
	
	save(list = ls(all=TRUE), file = "testWindows.RData")
}

scatter(WINDOW_SIZES, pPhaseWindow, "testWindows", xlab="Window size", ylab="Phaseability", main=params, type="b")
