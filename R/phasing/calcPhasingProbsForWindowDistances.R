calcPhasingProbsForWindowDistances <- function(distances, MAX_WINDOW_SIZE, FILE_NAME = NULL) {
	WINDOW_SIZES = 2:MAX_WINDOW_SIZE

	phaseProbsPositionWindow = matrix(data = NA, nrow=length(distances), ncol=length(WINDOW_SIZES))

	for (i in 1:length(distances)) {
		# Try to phase (i+1)-st position [relative to i] using varying window sizes:
		for (j in 1:length(WINDOW_SIZES)) {
			windowSize = WINDOW_SIZES[j]
			remainingSize = windowSize - 2  # exlcude i, i+1
			
			numOnLeft = i - 1
			numOnRight = (length(distances) + 1) - (i + 2) + 1
	
			if (numOnLeft <= numOnRight) {
				halfToUse = floor(remainingSize / 2) # skimp on the left [floor], and be generous with the right side
		                useOnLeft = min(halfToUse, numOnLeft)
		                useOnRight = min(remainingSize - useOnLeft, numOnRight)
			}
			else {
				halfToUse = ceiling(remainingSize / 2) # be generous with the right side [ceiling]
		                useOnRight = min(halfToUse, numOnRight)
		                useOnLeft = min(remainingSize - useOnRight, numOnLeft)
			}
			startInd = i - useOnLeft      # go left from position i
			stopInd = i + 1 + useOnRight  # go right from position i + 1
	
			usePositionRange = seq(from=startInd, to=stopInd, by=1)
			useDistancesRange = seq(from=startInd, to=stopInd-1, by=1)  # since there are N-1 distances between N consecutive positions
	
			phaseIndex = which(usePositionRange == i+1)
			if (length(phaseIndex) != 1) stop("NO phaseIndex!")
			windowDistances = distances[useDistancesRange]
	
			print(paste("Try to phase position ", i+1, " [relative to ", i, "] using positions: (", paste(usePositionRange, collapse=", "), "), windowDistances= (", paste(windowDistances, collapse=", "), "), [phaseIndex= ", phaseIndex, ", i=", i, "]", sep=""))
			p = pPhaseHetPairAtDistanceUsingDepthAndWindow(windowDistances, phaseIndex, meanDepth, nReadsToPhase, L, Im, Is)
			print(paste("phase prob: ", p, sep=""))
			phaseProbsPositionWindow[i, j] = p
		}
	
		if (!is.null(FILE_NAME)) {
	        	save(list = ls(all=TRUE), file = paste(FILE_NAME, ".RData", sep=""))
		}
	}
	
	phaseProbsPositionWindow
}
