L = 101

Fm = 392
Fs = 44

meanDepth = 65
nReadsToPhase = 1

params = paste("meanDepth= ", meanDepth, ", nReadsToPhase= ", nReadsToPhase, ", L= ", L, ", Fm= ", Fm, ", Fs= ", Fs, sep="")


DISTANCES = 0:1000
pPhaseHetPairAtDistWithRead = pDirectlyPhaseHetPairAtDistanceUsingDepth(meanDepth, nReadsToPhase, L, DISTANCES, Fm, Fs)

scatter(DISTANCES, pPhaseHetPairAtDistWithRead, "testSpecificDistances", xlab="Intra-het distance", ylab="Phaseability", main=params)



save(list = ls(all=TRUE), file = "testSpecificDistances.RData")
