L = 101

Fm = 392
Fs = 44

meanDepth = 65
nReadsToPhase = 1

params = paste("meanDepth= ", meanDepth, ", nReadsToPhase= ", nReadsToPhase, ", L= ", L, ", Fm= ", Fm, ", Fs= ", Fs, sep="")



MEAN_INTRA_HET_DISTANCES = seq(from=2, to=20002, by=50)
THETAS = meanIntraHetDistanceToTheta(MEAN_INTRA_HET_DISTANCES)
NUM_THETAS = length(THETAS)

pPhaseTheta = as.vector(matrix(data = -1, nrow = 1, ncol = NUM_THETAS))
for (i in 1:NUM_THETAS) {
	pPhaseTheta[i] = pDirectlyPhaseHetPair(meanDepth, nReadsToPhase, L, THETAS[i], Fm, Fs)
}
scatter(MEAN_INTRA_HET_DISTANCES, pPhaseTheta, "testIntraHetDistances", xlab="Mean intra-het distance", ylab="Phaseability", main=params, type="b")



save(list = ls(all=TRUE), file = "testIntraHetDistances.RData")
