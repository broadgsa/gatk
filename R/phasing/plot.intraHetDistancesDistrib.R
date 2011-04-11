theta = 10^-3
params = paste("theta= ", theta, sep="")

MIN_DIST = 1
MAX_DIST = 10^4
BY_DIST = 10
DISTANCES = seq(from=MIN_DIST, to=MAX_DIST+BY_DIST, by=BY_DIST)
freqAtLteDist = pHetPairLteDistance(DISTANCES, theta)

scatter(DISTANCES, freqAtLteDist, "intraHetDistancesDistrib", xlab="Intra-het distance", ylab="Cumulative Frequency", log="x", main=params)


save(list = ls(all=TRUE), file = "intraHetDistancesDistrib.RData")
