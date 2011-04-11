L = 101

Fm = 392
Fs = 44

meanDepth = 65
nReadsToPhase = 1

theta = 10^-3

params = paste("meanDepth= ", meanDepth, ", nReadsToPhase= ", nReadsToPhase, ", theta= ", theta, ", L= ", L, ", Fm= ", Fm, ", Fs= ", Fs, sep="")


MAX_NUM_DISTS = 10^4
distances = sampleIntraHetDistances(MAX_NUM_DISTS, theta)
print(paste("Using ", MAX_NUM_DISTS, " THEORETICAL distances...", sep=""))


MAX_WINDOW_SIZE = 10
FILE_NAME = "theoretical_window"

phaseWindowResult = calcPhasingProbsForWindowDistances(distances, MAX_WINDOW_SIZE, meanDepth, nReadsToPhase, L, Fm, Fs, FILE_NAME)
phaseProbsPositionWindow = phaseWindowResult$phaseProbsPositionWindow
WINDOW_SIZES = phaseWindowResult$WINDOW_SIZES

phaseProbsWindow = colMeans(phaseProbsPositionWindow)

scatter(WINDOW_SIZES, phaseProbsWindow, FILE_NAME, xlab="Window size", ylab="Mean theoretical phasing rate on empirical distances", main=params, type="b")
