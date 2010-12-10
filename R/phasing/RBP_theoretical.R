# For consecutive diploid het sites x and y, P(distance(x,y) = k)
# = P(site y is the first het site downstream of x at distance = k | het site x exists at its location).
# That is, het site x already "exists", and we want to know what the probability that the NEXT het site (y) is k bases away.
#
# pOneSiteIsHom = p(top chromosome is ref AND bottom chromosome is ref) + p(top chromosome is var AND bottom chromosome is var)
# = (1-theta)^2 + theta^2
#
# pOneSiteIsHet = p(top chromosome is ref AND bottom chromosome is var) + p(top chromosome is var AND bottom chromosome is ref)
# = (1-theta)*theta + theta*(1-theta) = 2*theta*(1-theta)
#
pHetPairAtDistance <- function(k, theta) {
	pOneSiteIsHet = 2 * theta * (1 - theta)
	dexp(k, pOneSiteIsHet)
}

# Since the geometric/exponential distribution is "memory-free", can simply multiply the (independent) probabilities for the distances:
pHetPairsAtDistances <- function(dists, theta) {
	prod(pHetPairAtDistance(dists, theta))
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

# For consecutive diploid het sites x and y, P(distance(x,y) <= k)
pHetPairLteDistance <- function(k, theta) {
	# Although the real minimum distance starts with 1 (geometric distribution), the exponential distribution approximation starts with 0:
	MIN_DISTANCE = 0

	Vectorize(function(maxDist) integrate(function(dist) pHetPairAtDistance(dist, theta), lower=MIN_DISTANCE, upper=maxDist)$value)(k)
}

# Probability (over locations of x on the read) that a paired-end read ALREADY covering site x [with 2 mates of length L and an insert size of i between them] will ALSO cover site y (k bases downstream of x):
#
# Assume that read is equally likely to cover x at any of the 2L positions, so uniform probability of 1/2L at each of them.
# P(read r covers (x,y) | r covers x, r = [L,i,L], distance(x,y) = k)
# = sum_p=1^p=L {1/2L * 1{k <= L-p OR L-p+i+1 <= k <= 2L+i-p}} + sum_p=1^p=L {1/2L * 1{k <= L-p}}
# = 1/2L * [2 * sum_p=1^p=L {1{k <= L-p}} + sum_p=1^p=L {1{L-p+i+1 <= k <= 2L+i-p}}]
# = 1/2L * [2 * max(0, L-k) + max(0, min(L, max(0, k-i)) - max(0, k-i-L))]
pReadWithSpecificInsertCanCoverHetPairAtDistance <- function(L, i, k) {
	pWithinSameMate = 2 * pmax(0, L - k)

	maxValueFor_p = pmin(L, pmax(0, k - i))
	minValueFor_p_minusOne = pmax(0, k - i - L)
	pInDifferentMates = pmax(0, maxValueFor_p - minValueFor_p_minusOne)

	(pWithinSameMate + pInDifferentMates) / (2*L)
}

# Probability of having an insert of size insertSize, where the insert sizes are normally distributed with mean Im and standard deviation Is:
pInsertSize <- function(insertSize, Im, Is) {
	dnorm(insertSize, mean = Im, sd = Is)
}

# Probability (over locations of x on the read, and insert sizes) that there could exist a paired-end read [with 2 mates of length L and an insert between them] covers both sites x and y (at distance k):
# Integral_from_0^to_INFINITY { pInsertSize(s, Im, Is) * pReadWithSpecificInsertCanCoverHetPairAtDistance(L, s, k) ds }
pReadCanCoverHetPairAtDistance <- function(L, k, Im, Is) {
	if (Is != 0) {
		pCoverageBySpecificInsert <- function(s) {pInsertSize(s, Im, Is) * pReadWithSpecificInsertCanCoverHetPairAtDistance(L, s, k)}

		MAX_NUM_SD = 10
		maxDistance = MAX_NUM_SD * Is
		minInsertSize = max(0, Im - maxDistance)
		maxInsertSize = Im + maxDistance
	
		integrate(pCoverageBySpecificInsert, lower=minInsertSize, upper=maxInsertSize)$value
	}
	else {# All reads have inserts of size exactly Im:
		pReadWithSpecificInsertCanCoverHetPairAtDistance(L, Im, k)
	}
}

# Probability (over locations of x on the read, insert sizes, and read depths) that there exist at least nReadsToPhase paired-end reads covering both sites x and y (at distance k):
# = Sum_from_d=0^to_d=2*meanDepth { p(having d reads | poisson with meanDepth) * p(there at least nReadsToPhase succeed in phasing x,y | given d reads in total) }
# p(having d reads | poisson with meanDepth) = dpois(d, meanDepth)
# p(there are at least nReadsToPhase that succeed in phasing x,y | given d reads in total) = pbinom(nReadsToPhase - 1, k, pReadCanCoverHetPairAtDistance(L, k, Im, Is), lower.tail = FALSE)
pDirectlyPhaseHetPairAtDistanceUsingDepth <- function(meanDepth, nReadsToPhase, L, k, Im, Is) {
	p = pReadCanCoverHetPairAtDistance(L, k, Im, Is)
	pAtLeastNreadsToPhaseGivenDepth <- function(d) pbinom(nReadsToPhase - 1, d, p, lower.tail = FALSE)
	pAtLeastNreadsToPhaseAndDepth	<- function(d) dpois(d, meanDepth) * pAtLeastNreadsToPhaseGivenDepth(d)
	
	minDepth = 0
	maxDepth = 2 * meanDepth
	sum(apply(as.matrix(minDepth:maxDepth), 1, pAtLeastNreadsToPhaseAndDepth))
}

pDirectlyPhaseHetPairAndDistanceUsingDepth <- function(meanDepth, nReadsToPhase, L, k, theta, Im, Is) {
	Vectorize(function(dist) pDirectlyPhaseHetPairAtDistanceUsingDepth(meanDepth, nReadsToPhase, L, dist, Im, Is) * pHetPairAtDistance(dist, theta))(k)
}

# Probability (over locations of x on the read, insert sizes, read depths, and het-het distances) that that there exist at least nReadsToPhase paired-end reads covering both sites x and y (where the distance between x and y is as per the geometric/exponential distribution):
pDirectlyPhaseHetPair <- function(meanDepth, nReadsToPhase, L, theta, Im, Is) {
	# Although the real minimum distance starts with 1 (geometric distribution), the exponential distribution approximation starts with 0:
	MIN_DISTANCE = 0
	MAX_DISTANCE = Inf

	integrate(function(k) pDirectlyPhaseHetPairAndDistanceUsingDepth(meanDepth, nReadsToPhase, L, k, theta, Im, Is), lower=MIN_DISTANCE, upper=MAX_DISTANCE, subdivisions=1000)$value
}

# Probability (over locations of sites on reads, insert sizes, and read depths) that paired-end reads can TRANSITIVELY phase phaseIndex relative to phaseIndex - 1, given a window of length(windowDistances)+1 het sites at distances given by windowDistances (where an edge in the transitive path requires at least nReadsToPhase reads):
pPhaseHetPairAtDistanceUsingDepthAndWindow <- function(windowDistances, phaseIndex, meanDepth, nReadsToPhase, L, Im, Is) {
	n = length(windowDistances) + 1  # the window size
	if (phaseIndex < 2 || phaseIndex > n) {
		stop("phaseIndex < 2 || phaseIndex > n")
	}

	# A. Pre-compute the upper diagonal of square matrix of n CHOOSE 2 values of:
	# pDirectlyPhaseHetPairAtDistanceUsingDepth(meanDepth, nReadsToPhase, L, dist(i,j), Im, Is)
	#
	# NOTE that the probabilities of phasing different pairs are NOT truly independent, but assume this for convenience...
	#
	pPhasePair = matrix(data = 0, nrow = n, ncol = n)
	for (i in seq(from=1, to=n-1, by=1)) {
		for (j in seq(from=i+1, to=n, by=1)) {
			dist = distanceBetweenPair(i, j, windowDistances)
			#print(paste("distanceBetweenPair(", i, ", ", j, ", windowDistances) = ", dist, sep=""))
			
			pPhaseIandJ = pDirectlyPhaseHetPairAtDistanceUsingDepth(meanDepth, nReadsToPhase, L, dist, Im, Is)
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

pPhaseHetPairAndDistancesUsingDepthAndWindow <- function(windowDistances, phaseIndex, meanDepth, nReadsToPhase, L, Im, Is, theta) {
	p = pPhaseHetPairAtDistanceUsingDepthAndWindow(windowDistances, phaseIndex, meanDepth, nReadsToPhase, L, Im, Is) * pHetPairsAtDistances(windowDistances, theta)

	#print(paste(p, " = pPhaseHetPairAndDistancesUsingDepthAndWindow(windowDistances= (", paste(windowDistances, collapse=", "), "), phaseIndex= ", phaseIndex, ", meanDepth= ", meanDepth, ", nReadsToPhase= ", nReadsToPhase, ", L= ", L, ", Im= ", Im, ", Is= ", Is, ", theta= ", theta, ") * pHetPairsAtDistances(windowDistances= ", paste(windowDistances, collapse=", "), ", theta= ", theta, ")", sep=""))

	p
}

# Probability (over locations of sites on reads, insert sizes, and read depths) that paired-end reads can TRANSITIVELY phase phaseIndex relative to phaseIndex - 1, given a window of n het sites at distances distributed as determined by theta (where an edge in the transitive path requires at least nReadsToPhase reads):
pDirectlyPhaseHetPairUsingWindow <- function(meanDepth, nReadsToPhase, L, theta, Im, Is, n, phaseIndex) {
	if (n < 2) {
		stop("n < 2")
	}
	ndim = n-1

	integrandFunction <- function(windowDistances) {pPhaseHetPairAndDistancesUsingDepthAndWindow(windowDistances, phaseIndex, meanDepth, nReadsToPhase, L, Im, Is, theta)}

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
