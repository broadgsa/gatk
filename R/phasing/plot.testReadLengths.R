theta = 10^-3

Fm_BASE = 392 - 2 * 101  # The mean insert size == 190
Fs = 44

meanDepth = 65
nReadsToPhase = 1

params = paste("meanDepth= ", meanDepth, ", nReadsToPhase= ", nReadsToPhase, ", theta= ", theta, ", Fm_BASE= ", Fm_BASE, ", Fs= ", Fs, sep="")


READ_LENGTHS = seq(from=30, to=1000, by=10)
NUM_READ_LENGTHS = length(READ_LENGTHS)

pPhaseReadLength = as.vector(matrix(data = -1, nrow = 1, ncol = NUM_READ_LENGTHS))
for (i in 1:NUM_READ_LENGTHS) {
	Fm = Fm_BASE + 2 * READ_LENGTHS[i]
	pPhaseReadLength[i] = pDirectlyPhaseHetPair(meanDepth, nReadsToPhase, READ_LENGTHS[i], theta, Fm, Fs)
}
scatter(READ_LENGTHS, pPhaseReadLength, "testReadLengths", xlab="Read length", ylab="Phaseability", main=params, type="b")



save(list = ls(all=TRUE), file = "testReadLengths.RData")
