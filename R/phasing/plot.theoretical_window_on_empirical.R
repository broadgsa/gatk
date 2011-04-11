L = 101

Fm = 392
Fs = 44

meanDepth = 65
nReadsToPhase = 1

params = paste("meanDepth= ", meanDepth, ", nReadsToPhase= ", nReadsToPhase, ", L= ", L, ", Fm= ", Fm, ", Fs= ", Fs, sep="")


distances = scan("~fromer/storage/phase.NA12878/COMPLETE_LIST.het_distances.txt", what=list(dist=0))
distances = distances$dist

MAX_NUM_DISTS = 10^4
NUM_DISTS_TO_USE = min(MAX_NUM_DISTS, length(distances))
distances = distances[1:NUM_DISTS_TO_USE]
print(paste("Using ", NUM_DISTS_TO_USE, " EMPIRICAL distances...", sep=""))


MAX_WINDOW_SIZE = 10
FILE_NAME = "theoretical_window_on_empirical"

phaseWindowResult = calcPhasingProbsForWindowDistances(distances, MAX_WINDOW_SIZE, meanDepth, nReadsToPhase, L, Fm, Fs, FILE_NAME)
phaseProbsPositionWindow = phaseWindowResult$phaseProbsPositionWindow
WINDOW_SIZES = phaseWindowResult$WINDOW_SIZES

phaseProbsWindow = colMeans(phaseProbsPositionWindow)

scatter(WINDOW_SIZES, phaseProbsWindow, FILE_NAME, xlab="Window size", ylab="Mean theoretical phasing rate on empirical distances", main=params, type="b")
