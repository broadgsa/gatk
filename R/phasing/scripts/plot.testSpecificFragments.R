L = 76
k = 75
params = paste("L= ", L, ", k= ", k, sep="")

FRAGMENT_SIZES = 0:100 + 2 * L
pCoverHetPairWithRead = pPairedEndReadsOfSpecificFragmentCanCoverHetPairAtDistance(L, FRAGMENT_SIZES, k)

scatter(FRAGMENT_SIZES, pCoverHetPairWithRead, "testSpecificFragments", xlab="Fragment size", ylab="Probability of covering het pair", main=params)
