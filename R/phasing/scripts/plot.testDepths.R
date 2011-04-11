theta = 10^-3

Fm_BASE = 392 - 2 * 101  # The mean insert size == 190
Fs = 44

nReadsToPhase = 1

params = paste("nReadsToPhase= ", nReadsToPhase, ", theta= ", theta, ", Fm_BASE= ", Fm_BASE, ", Fs= ", Fs, sep="")



MEAN_DEPTHS = 0:65
NUM_DEPTHS = length(MEAN_DEPTHS)

READ_LENGTHS = c(18, 36, 76, 101, 125, 150, 175, 200, 400, 800, 1000)
READ_LENGTHS = rev(READ_LENGTHS)
NUM_READ_LENGTHS = length(READ_LENGTHS)

depthsX = list()
depthsY = list()
depthsLeg = vector()

for (i in 1:NUM_READ_LENGTHS) {
	pPhaseDepth = as.vector(matrix(data = -1, nrow = 1, ncol = NUM_DEPTHS))
	Fm = Fm_BASE + 2 * READ_LENGTHS[i]
	for (j in 1:NUM_DEPTHS) {
		pPhaseDepth[j] = pDirectlyPhaseHetPair(MEAN_DEPTHS[j], nReadsToPhase, READ_LENGTHS[i], theta, Fm, Fs)
	}
	depthsX[i] = list(MEAN_DEPTHS)
	depthsY[i] = list(pPhaseDepth)
	depthsLeg[i] = paste("L= ", READ_LENGTHS[i], sep="")
}

scatter(depthsX, depthsY, "testDepths", xlab="Mean depth", ylab="Phaseability", main=params, leg=depthsLeg, legPos="topleft", width=14, height=7, type="b")




save(list = ls(all=TRUE), file = "testDepths.RData")
