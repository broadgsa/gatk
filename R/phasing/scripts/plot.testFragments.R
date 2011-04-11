theta = 10^-3

L = 101

meanDepth = 65
nReadsToPhase = 1

params = paste("meanDepth= ", meanDepth, ", nReadsToPhase= ", nReadsToPhase, ", L= ", L, ", theta= ", theta, sep="")



MEAN_SIZES = seq(1,2000,20)
STD_SIZES = seq(0,200,5)


testFragments = matrix(nrow=length(MEAN_SIZES), ncol=length(STD_SIZES))
for (i in 1:length(MEAN_SIZES)) {
	test_mean_fragment_size = MEAN_SIZES[i]
	print(paste("test_mean_fragment_size: ", test_mean_fragment_size, sep=""))
	for (j in 1:length(STD_SIZES)) {
		test_std_fragment_size = STD_SIZES[j]
		print(paste("test_std_fragment_size: ", test_std_fragment_size, sep=""))

		testFragments[i,j] = pDirectlyPhaseHetPair(meanDepth, nReadsToPhase, L, theta, test_mean_fragment_size, test_std_fragment_size)
	}
}


pdf('testFragments.pdf')

library(gplots)
heatmap.2(testFragments, ylab = "Mean fragment size", xlab = "Standard deviation fragment size", labRow = MEAN_SIZES, labCol = STD_SIZES, Rowv = NA, Colv = NA, dendrogram = "none", scale="none", revC = FALSE, density.info="none", trace="none", main=params)

library(scatterplot3d)
xMeans = as.vector(t(matrix(rep.int(MEAN_SIZES, length(STD_SIZES)), ncol = length(STD_SIZES))))
yStds = rep.int(STD_SIZES, length(MEAN_SIZES))
zPhaseRate = as.vector(t(testFragments))
scatterplot3d(xMeans, yStds, zPhaseRate, xlab = "Mean fragment size", ylab = "Standard deviation fragment size", zlab = "Phasing rate", main=params)

bestCombo = which.max(zPhaseRate)
print(paste("For ", params, ", BEST choice gives phaseability of ", zPhaseRate[bestCombo], " using mean fragment = ", xMeans[bestCombo], ", std. fragment = ", yStds[bestCombo], sep = ""))
dev.off()




save(list = ls(all=TRUE), file = "testFragments.RData")
