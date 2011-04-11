#
#options(warn=2)
#options(error=recover)
#

HALF = high_dimensional_integrate(1, -200, 0, dnorm)
print(paste("Should be ~ HALF: ", HALF, sep=""))


k = 75
#theta = 10^-2
theta = 10^-3

p = pHetPairLteDistance(k, theta)
print(paste(p, " = pHetPairLteDistance(k= ", k, ", theta= ", theta, ")", sep=""))


L = 76
fragmentSize = 452


p = pPairedEndReadsOfSpecificFragmentCanCoverHetPairAtDistance(L, fragmentSize, k)
print(paste(p, " = pPairedEndReadsOfSpecificFragmentCanCoverHetPairAtDistance(L= ", L, ", fragmentSize= ", fragmentSize, ", k= ", k, ")", sep=""))

Fm = 392
Fs = 44

p = pFragmentSize(300, Fm, Fs)
print(paste(p, " = pFragmentSize(300, Fm= ", Fm, ", Fs= ", Fs, ")", sep=""))


p = pFragmentsReadsCanCoverHetPairAtDistance(L, k, Fm, Fs)
print(paste(p, " = pFragmentsReadsCanCoverHetPairAtDistance(L= ", L, ", k= ", k, ", Fm= ", Fm, ", Fs= ", Fs, ")", sep=""))


meanDepth = 65
nReadsToPhase = 1
p = pDirectlyPhaseHetPairAtDistanceUsingDepth(meanDepth, nReadsToPhase, L, k, Fm, Fs)
print(paste(p, " = pDirectlyPhaseHetPairAtDistanceUsingDepth(meanDepth= ", meanDepth, ", nReadsToPhase= ", nReadsToPhase, ", L= ", L, ", k= ", k, ", Fm= ", Fm, ", Fs= ", Fs, ")", sep=""))


p = pDirectlyPhaseHetPair(meanDepth, nReadsToPhase, L, theta, Fm, Fs)
print(paste(p, " = pDirectlyPhaseHetPair(meanDepth= ", meanDepth, ", nReadsToPhase= ", nReadsToPhase, ", L= ", L, ", theta= ", theta, ", Fm= ", Fm, ", Fs= ", Fs, ")", sep=""))


windowDistances = c(100, 100, 100, 100, 100)
phaseIndex = 2
p = pPhaseHetPairAtDistanceUsingDepthAndWindow(windowDistances, phaseIndex, meanDepth, nReadsToPhase, L, Fm, Fs)
print(paste(p, " = pPhaseHetPairAtDistanceUsingDepthAndWindow(windowDistances= (", paste(windowDistances, collapse=", "), "), phaseIndex= ", phaseIndex, ", meanDepth= ", meanDepth, ", nReadsToPhase= ", nReadsToPhase, ", L= ", L, ", Fm= ", Fm, ", Fs= ", Fs, ")", sep=""))



traceback()
warnings()
