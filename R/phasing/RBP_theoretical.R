# pOneSiteIsHom = p(top chromosome is ref AND bottom chromosome is ref) + p(top chromosome is var AND bottom chromosome is var)
# = (1-theta)^2 + theta^2
#
# pOneSiteIsHet = p(top chromosome is ref AND bottom chromosome is var) + p(top chromosome is var AND bottom chromosome is ref)
# = (1-theta)*theta + theta*(1-theta) = 2*theta*(1-theta)
pOneSiteIsHet <- function(theta) {
	2 * theta * (1 - theta)
}

# p = 2 * theta * (1 - theta)
# and mean intra-het distance = 1/p, or d = 1/p
# or: p = 1/d
# or: 2 * theta * (1 - theta) = 1/d
# theta * (1 - theta) = 1/2d
# - theta^2 + theta - 1/2d = 0
#
# Using the quadratic equation:
# (- b + (b^2 - 4*a*c)^0.5) / 2a
# (-1 + (1 - 2/d)^0.5) / -2
meanIntraHetDistanceToTheta <- function(d) {
	(-1 + (1 - 2/d)^0.5) / -2
}

# For consecutive diploid het sites x and y, P(distance(x,y) = k)
# = P(site y is the first het site downstream of x at distance = k | het site x exists at its location).
# That is, het site x already "exists", and we want to know what the probability that the NEXT het site (y) is k bases away.
pHetPairAtDistance <- function(k, theta) {
	pOneSiteIsHetTheta = pOneSiteIsHet(theta)
	dexp(k, pOneSiteIsHetTheta)
}

# Since the geometric/exponential distribution is "memory-free", can simply multiply the (independent) probabilities for the distances:
pHetPairsAtDistances <- function(dists, theta) {
	prod(pHetPairAtDistance(dists, theta))
}

# Sample numDists distances from the intra-het distance distribution.
# [since the geometric/exponential distribution is "memory-free", can simply **independently** sample from the distribution]:
sampleIntraHetDistances <- function(numDists, theta) {
	pOneSiteIsHetTheta = pOneSiteIsHet(theta)
	ceiling(rexp(numDists, pOneSiteIsHetTheta))  # round up to get whole-number distances starting from 1
}

# For consecutive diploid het sites x and y, P(distance(x,y) <= k)
pHetPairLteDistance <- function(k, theta) {
	# Although the real minimum distance starts with 1 (geometric distribution), the exponential distribution approximation starts with 0:
	MIN_DISTANCE = 0

	Vectorize(function(maxDist) integrate(function(dist) pHetPairAtDistance(dist, theta), lower=MIN_DISTANCE, upper=maxDist)$value)(k)
}

# Probability (over locations of x on the read) that a paired-end read ALREADY covering site x [with 2 mates of length L reading a fragment of length F] will ALSO cover site y (k bases downstream of x):
#
# If read 1 in mate spans [s1, e1] and read 2 spans [s2, e2], where length(read 1) = e1 - s1 + 1 = length(read 2) = e2 - s2 + 1 = L, then i = s2 - e1 - 1 [BY DEFINITION of i].
# i == "insert size" is DEFINED AS: F - 2 * L
#
#
# FOR i >= 0:
#
# Assume that read is equally likely to cover x at any of the 2L positions, so uniform probability of 1/2L at each of them.
# P(read r covers (x,y) | r covers x, r = [L,i,L], distance(x,y) = k)
# = sum_p=1^p=L {1/2L * 1{k <= L-p OR L-p+i+1 <= k <= 2L+i-p}} + sum_p=1^p=L {1/2L * 1{k <= L-p}}
# = 1/2L * [2 * sum_p=1^p=L {1{k <= L-p}} + sum_p=1^p=L {1{L-p+i+1 <= k <= 2L+i-p}}]
# = 1/2L * [2 * max(0, L-k) + max(0, min(L, max(0, k-i)) - max(0, k-i-L))]
#
#
pPairedEndReadsOfSpecificFragmentCanCoverHetPairAtDistance <- function(L, F, k) {
	if (min(F) < 1) {
		stop("Cannot have fragments of size < 1")
	}

	# if F < L, then set the effective read length to be F:
	L = pmin(L, F)

	i = F - 2 * L
	#print(paste("pPairedEndReadsOfSpecificFragmentCanCoverHetPairAtDistance(L= (", paste(L, collapse=", "), "), F= (", paste(F, collapse=", "), "), k= (", paste(k, collapse=", "), ")), i= (", paste(i, collapse=", "), ")", sep=""))

	# If i < 0, then ASSUMING that overlapping region is identical, we can "pretend" to have 2 reads of length L and L+i, with no insert between them.
	# Otherwise, leave i alone and L1 = L2 = L:
	L1 = L
	L2 = L + pmin(0, i)   # set effective length of second read to L+i if i < 0
	i = pmax(0, i) # set effective insert size to be >= 0


	pWithinSameMate = pmax(0, L1 - k) + pmax(0, L2 - k)

	#maxValueFor_p = pmin(L1, pmax(0, k - i))
	#minValueFor_p_minusOne = pmax(0, k - i - L2)

	maxValueFor_p = pmin(L1, L1 + L2 + i - k)
	minValueFor_p_minusOne = pmax(0, L1 - k + i)
	pInDifferentMates = pmax(0, maxValueFor_p - minValueFor_p_minusOne)

	(pWithinSameMate + pInDifferentMates) / (L1 + L2)
}

# Probability of having a fragment of size fragmentSize, where the fragment sizes are normally distributed with mean Fm and standard deviation Fs:
pFragmentSize <- function(fragmentSize, Fm, Fs) {
	dnorm(fragmentSize, mean = Fm, sd = Fs)
}

# Probability (over locations of x on the read, and fragment sizes) that there could exist a paired-end read [with 2 mates of length L covering a fragment] covers both sites x and y (at distance k):
# Integral_from_0^to_INFINITY { pFragmentSize(s, Fm, Fs) * pPairedEndReadsOfSpecificFragmentCanCoverHetPairAtDistance(L, s, k) ds }
pFragmentsReadsCanCoverHetPairAtDistance <- function(L, k, Fm, Fs) {
	if (Fs != 0) {
		pCoverageBySpecificFragment <- function(s) {pFragmentSize(s, Fm, Fs) * pPairedEndReadsOfSpecificFragmentCanCoverHetPairAtDistance(L, s, k)}

		MAX_NUM_SD = 10
		maxDistance = MAX_NUM_SD * Fs
		minFragmentSize = max(1, Fm - maxDistance) # NOT meaningful to have fragment size < 1
		maxFragmentSize = Fm + maxDistance
	
		integrate(pCoverageBySpecificFragment, lower=minFragmentSize, upper=maxFragmentSize)$value
	}
	else {# All fragments are of size exactly Fm:
		pPairedEndReadsOfSpecificFragmentCanCoverHetPairAtDistance(L, Fm, k)
	}
}

# Probability (over locations of x on the read, fragment sizes, and read depths) that there exist at least nReadsToPhase paired-end reads covering both sites x and y (at distance k):
# = Sum_from_d=0^to_d=2*meanDepth { p(having d reads | poisson with meanDepth) * p(there at least nReadsToPhase succeed in phasing x,y | given d reads in total) }
# p(having d reads | poisson with meanDepth) = dpois(d, meanDepth)
# p(there are at least nReadsToPhase that succeed in phasing x,y | given d reads in total) = pbinom(nReadsToPhase - 1, k, pFragmentsReadsCanCoverHetPairAtDistance(L, k, Fm, Fs), lower.tail = FALSE)
pDirectlyPhaseHetPairAtDistanceUsingDepth_SINGLE_k <- function(meanDepth, nReadsToPhase, L, k, Fm, Fs) {
	THRESH = 10^-8
	p = pFragmentsReadsCanCoverHetPairAtDistance(L, k, Fm, Fs)

	# deal with numerical issues:
	if (abs(1 - p) < THRESH) {
		p = 1
	}
	else if (abs(p) < THRESH) {
		p = 0
	}

	pAtLeastNreadsToPhaseGivenDepth <- function(d) pbinom(nReadsToPhase - 1, d, p, lower.tail = FALSE)
	pAtLeastNreadsToPhaseAndDepth	<- function(d) dpois(d, meanDepth) * pAtLeastNreadsToPhaseGivenDepth(d)
	
	minDepth = 0
	maxDepth = 2 * meanDepth
	sum(apply(as.matrix(minDepth:maxDepth), 1, pAtLeastNreadsToPhaseAndDepth))
}

pDirectlyPhaseHetPairAtDistanceUsingDepth <- function(meanDepth, nReadsToPhase, L, k, Fm, Fs) {
	Vectorize(function(dist) pDirectlyPhaseHetPairAtDistanceUsingDepth_SINGLE_k(meanDepth, nReadsToPhase, L, dist, Fm, Fs))(k)
}

pDirectlyPhaseHetPairAndDistanceUsingDepth <- function(meanDepth, nReadsToPhase, L, k, theta, Fm, Fs) {
	pDirectlyPhaseHetPairAtDistanceUsingDepth(meanDepth, nReadsToPhase, L, k, Fm, Fs) * pHetPairAtDistance(k, theta)
}

# Probability (over locations of x on the read, fragment sizes, read depths, and het-het distances) that that there exist at least nReadsToPhase paired-end reads covering both sites x and y (where the distance between x and y is as per the geometric/exponential distribution):
pDirectlyPhaseHetPair <- function(meanDepth, nReadsToPhase, L, theta, Fm, Fs) {
	# Although the real minimum distance starts with 1 (geometric distribution), the exponential distribution approximation starts with 0:
	MIN_DISTANCE = 0
	MAX_DISTANCE = Inf

	iRes = integrate(function(k) pDirectlyPhaseHetPairAndDistanceUsingDepth(meanDepth, nReadsToPhase, L, k, theta, Fm, Fs), lower=MIN_DISTANCE, upper=MAX_DISTANCE, subdivisions=1000, stop.on.error = FALSE)
	if (iRes$message != "OK") {
		print(paste("DISTANCE INTEGRATION WARNING: ", iRes$message, sep=""))
	}
	iRes$value
}

# Probability (over locations of sites on reads, fragment sizes, and read depths) that paired-end reads can TRANSITIVELY phase phaseIndex relative to phaseIndex - 1, given a window of length(windowDistances)+1 het sites at distances given by windowDistances (where an edge in the transitive path requires at least nReadsToPhase reads):
pPhaseHetPairAtDistanceUsingDepthAndWindow <- function(windowDistances, phaseIndex, meanDepth, nReadsToPhase, L, Fm, Fs, MIN_PATH_PROB = 10^-6) {
	n = length(windowDistances) + 1  # the window size
	if (phaseIndex < 2 || phaseIndex > n) {
		stop("phaseIndex < 2 || phaseIndex > n")
	}
	#print(paste("windowDistances= (", paste(windowDistances, collapse=", "), ")", sep=""))

	# A. Pre-compute the upper diagonal of square matrix of n CHOOSE 2 values of:
	# pDirectlyPhaseHetPairAtDistanceUsingDepth(meanDepth, nReadsToPhase, L, dist(i,j), Fm, Fs)
	#
	# NOTE that the probabilities of phasing different pairs are NOT truly independent, but assume this for convenience...
	#
	pPhasePair = matrix(data = 0, nrow = n, ncol = n)
	for (i in seq(from=1, to=n-1, by=1)) {
		for (j in seq(from=i+1, to=n, by=1)) {
			dist = distanceBetweenPair(i, j, windowDistances)
			#print(paste("distanceBetweenPair(", i, ", ", j, ", windowDistances) = ", dist, sep=""))
			
			pPhaseIandJ = pDirectlyPhaseHetPairAtDistanceUsingDepth(meanDepth, nReadsToPhase, L, dist, Fm, Fs)
			pPhasePair[i, j] = pPhaseIandJ
			pPhasePair[j, i] = pPhaseIandJ
		}
	}
	#print(pPhasePair)

	# B. We need to consider ALL possible paths from phaseIndex - 1 ---> phaseIndex
	# There are:  sum_i=0^to_n-2  {n-2 CHOOSE i * i!}   such paths.
	# Multiply the phasing probs along the path, and sum over all such paths:
	#
	startNode = phaseIndex - 1
	endNode = phaseIndex

	possibleIntermediateNodes = vector()
	if (startNode > 1) possibleIntermediateNodes = c(possibleIntermediateNodes, seq(from=1, to=startNode-1, by=1))
	if (endNode < n) possibleIntermediateNodes = c(possibleIntermediateNodes, seq(from=endNode+1, to=n, by=1))
	#print(paste("possibleIntermediateNodes= {", paste(possibleIntermediateNodes, collapse=", "), "}", sep=""))

	pWindowNotPhasing = 1
	library(gtools)
	for (subset in powerSet(length(possibleIntermediateNodes))) {
		subset = possibleIntermediateNodes[subset]
		#print((paste("subset = {", paste(subset, collapse=", "), "}", sep="")))

		if (length(subset) == 0) {
			paths = c()
		}
		else {
			paths = permutations(length(subset), length(subset), v=subset)
		}
		# Add on the start and the end:
		paths = cbind(startNode, paths, endNode)

		for (i in 1:nrow(paths)) {
			path = paths[i,]
			pSpecificPathPhases = 1
			for (j in seq(from=1, to=length(path)-1, by=1)) {
				pSpecificPathPhases = pSpecificPathPhases * pPhasePair[path[j], path[j+1]]
				if (pSpecificPathPhases < MIN_PATH_PROB) { # Do a "bounded" calculation [any path that is ALREADY of low probability can be discarded]:
					#print(paste("pSpecificPathPhases= ", pSpecificPathPhases, sep=""))
					pSpecificPathPhases = 0
					break
				}
			}
			pWindowNotPhasing = pWindowNotPhasing * (1 - pSpecificPathPhases)

			#print((paste("path = (", paste(path, collapse=", "), "), pSpecificPathPhases= ", pSpecificPathPhases, sep="")))
		}
	}

	1 - pWindowNotPhasing
}

# distance(i,j) = distance(i,i+1) + ... + distance(j-1,j), where distance(i,i+1) is given by windowDistances(i):
distanceBetweenPair <- function(i, j, windowDistances) {
	if (i > j) {
		tmp = i
		i = j
		j = tmp
	}
	if (i < 1 || j > length(windowDistances) + 1) {
		stop(paste(i, " = i < 1 || ", j, " = j > length(windowDistances) + 1 = ", length(windowDistances) + 1, sep=""))
	}
	
	sum(windowDistances[i:(j-1)])
}

# n = size of set for which power set is to be returned
powerSet <- function(n) {
	library(sfsmisc)

	subsets = list()
	for (i in seq(from=0, to=(2^n)-1, by=1)) {
		subsets[i+1] = list(which(digitsBase(i, base = 2, ndigits = n) == 1))
	}
	subsets
}

pPhaseHetPairAndDistancesUsingDepthAndWindow <- function(windowDistances, phaseIndex, meanDepth, nReadsToPhase, L, Fm, Fs, theta) {
	p = pPhaseHetPairAtDistanceUsingDepthAndWindow(windowDistances, phaseIndex, meanDepth, nReadsToPhase, L, Fm, Fs) * pHetPairsAtDistances(windowDistances, theta)

	#print(paste(p, " = pPhaseHetPairAndDistancesUsingDepthAndWindow(windowDistances= (", paste(windowDistances, collapse=", "), "), phaseIndex= ", phaseIndex, ", meanDepth= ", meanDepth, ", nReadsToPhase= ", nReadsToPhase, ", L= ", L, ", Fm= ", Fm, ", Fs= ", Fs, ", theta= ", theta, ") * pHetPairsAtDistances(windowDistances= ", paste(windowDistances, collapse=", "), ", theta= ", theta, ")", sep=""))

	p
}

# Probability (over locations of sites on reads, fragment sizes, and read depths) that paired-end reads can TRANSITIVELY phase phaseIndex relative to phaseIndex - 1, given a window of n het sites at distances distributed as determined by theta (where an edge in the transitive path requires at least nReadsToPhase reads):
pDirectlyPhaseHetPairUsingWindow <- function(meanDepth, nReadsToPhase, L, theta, Fm, Fs, n, phaseIndex) {
	if (n < 2) {
		stop("n < 2")
	}
	ndim = n-1

	integrandFunction <- function(windowDistances) {pPhaseHetPairAndDistancesUsingDepthAndWindow(windowDistances, phaseIndex, meanDepth, nReadsToPhase, L, Fm, Fs, theta)}

	MIN_DISTANCE = 0

	#
	#MAX_DISTANCE = Inf
	#
	MAX_TAIL_PROB = 10^-6
	MAX_DISTANCE = 7500  # Only 3e-07 [=  1 - pHetPairLteDistance(7500, 10^-3)] of the het-het pairs are at a distance > 7500
	while (1 - pHetPairLteDistance(MAX_DISTANCE, theta) > MAX_TAIL_PROB) {
		MAX_DISTANCE = MAX_DISTANCE * 2
	}

	lower = as.vector(matrix(data=MIN_DISTANCE, nrow=1, ncol=ndim))
	upper = as.vector(matrix(data=MAX_DISTANCE, nrow=1, ncol=ndim))

	N = 10^4 * ndim^2
	high_dimensional_integrate(ndim, lower, upper, integrandFunction, N, DEBUG = TRUE, PRINT_EVERY = 10^2)
}

# Use the simplest version of the Monte Carlo method to integrate over a high-dimensional function:
high_dimensional_integrate <- function(ndim, lower, upper, integrandFunction, N = 10^4, DEBUG = FALSE, PRINT_EVERY = 10^3) {
	rectangularVolume = prod(upper - lower)

	sum = 0
	for (i in 1:N) {
		randVals = as.vector(matrix(data = NA, nrow=1, ncol=ndim))
		for (j in 1:ndim) {
			randVals[j] = runif(1, min=lower[j], max=upper[j])
		}
		#print(randVals)

		evalFuncVal = integrandFunction(randVals)
		sum = sum + evalFuncVal

		if (DEBUG && (i-1) %% PRINT_EVERY == 0) {
			estimate = rectangularVolume * (sum / i)
			print(paste("high_dimensional_integrate: iteration ", i, ", estimate= ", estimate, sep=""))
		}
	}
	rectangularVolume * (sum / N)
}

middleOfWindowIndex <- function(windowSize) {
	floor(windowSize/2 + 1)
}
