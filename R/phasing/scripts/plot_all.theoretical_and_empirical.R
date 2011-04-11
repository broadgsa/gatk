NUM_chr1_HET_SITES = as.integer(system("grep -c 'chr1:' ~fromer/storage/phase.NA12878/COMPLETE_LIST.het_sites.interval_list", intern=TRUE))
NUM_chr1_PHASEABLE_HET_SITES = NUM_chr1_HET_SITES - 1  # since can't phase the first het site


#
#USE_EMPIRICAL_WINDOWS = c(10, 2)
#
USE_EMPIRICAL_WINDOWS = c(2)


TWO_COLORS = c("red", "darkgreen")


######################################################################
# Phasing as a function of SPECIFIC intra-het distances:
######################################################################
load("testSpecificDistances.RData")

MAX_DISTANCE = 10^3
PQ_PHASING_THRESH = 10.0

distances = list()
phaseRateDistances = list()
distancesLeg = vector()

for (nextIndex in 1:length(USE_EMPIRICAL_WINDOWS)) {
	n = USE_EMPIRICAL_WINDOWS[nextIndex]
	n_locDistancePQReadsWindow <- scan(paste("~fromer/storage/phase.NA12878/phase_all_chr.n_", n, ".NA12878", ".locus_distance_PQ_numReads_windowSize.txt", sep=""), what=list(loci="", distance=0, PQ=0, reads=0, window=0))
	n_distance <- n_locDistancePQReadsWindow$distance
	n_PQ <- n_locDistancePQReadsWindow$PQ

	distanceVector = sort(unique(n_distance))
	distanceVector = distanceVector[which(distanceVector <= MAX_DISTANCE)]
	numDists = length(distanceVector)

	phasedFractionVector = as.vector(matrix(data=-1, nrow=1, ncol=numDists))
	
	print(paste("numDists= ", numDists, sep=""))
	print(paste(distanceVector, collapse=", "))

	for (i in 1:numDists) {
		d = distanceVector[i]
		print(paste("d= ", d, sep=""))

		dInds = which(n_distance == d)
		phasedFractionVector[i] = length(which(n_PQ[dInds] >= PQ_PHASING_THRESH)) / length(dInds)
	}

	distances[nextIndex] = list(distanceVector)
	phaseRateDistances[nextIndex] = list(phasedFractionVector)
	distancesLeg[nextIndex] = paste("HiSeq (window = ", n, ")", sep="")
}

nextIndex = nextIndex+1
distances[nextIndex] = list(DISTANCES)
phaseRateDistances[nextIndex] = list(pPhaseHetPairAtDistWithRead)
distancesLeg[nextIndex] = "Theoretical (window = 2)" # params

scatter(distances, phaseRateDistances, "specific_distances.theoretical_empirical", xlab="Intra-het distance", ylab="Phaseability", leg=distancesLeg, legPos="topright", width=14, height=7, type="b", col=TWO_COLORS)



######################################################################
# Phasing as a function of depth:
######################################################################
load("testDepths.RData")

depths = list()
phaseRateDepths = list()
depthsLeg = vector()

for (nextIndex in 1:length(USE_EMPIRICAL_WINDOWS)) {
	n = USE_EMPIRICAL_WINDOWS[nextIndex]
	RGdocPhasedConsistentSwitch = scan(paste("~fromer/storage/downsampled_phasing.NA12878.HiSeq/RG.DoC_phased_consistent_switch.chr1.n_", n, ".txt", sep=""), what=list(RGdoc=0, phased=0, consistentPhased=0, switch=0.0))
	depths[nextIndex] = list(RGdocPhasedConsistentSwitch$RGdoc)
	phaseRateDepths[nextIndex] = list(RGdocPhasedConsistentSwitch$phased / NUM_chr1_PHASEABLE_HET_SITES)
	depthsLeg[nextIndex] = paste("Down-sampled HiSeq (window = ", n, ")", sep="")
}

nextIndex = nextIndex+1
useLength = which(READ_LENGTHS == 101)
depths[nextIndex] = depthsX[useLength]
phaseRateDepths[nextIndex] = depthsY[useLength]
depthsLeg[nextIndex] = "Theoretical (window = 2)" # params

scatter(depths, phaseRateDepths, "depths.theoretical_empirical", xlab="Mean depth", ylab="Phaseability", leg=depthsLeg, legPos="topleft", width=14, height=7, type="b", col=TWO_COLORS)



######################################################################
# Distribution of intra-het distances:
######################################################################
load("intraHetDistancesDistrib.RData")

empiricalIntraHetDistances = read.table("~fromer/storage/phase.NA12878/COMPLETE_LIST.het_distances.txt")$V1
empiricalIntraHetDistances[which(empiricalIntraHetDistances >= MAX_DIST)] = MAX_DIST

empiricalIntraHetDistancesHist = hist(empiricalIntraHetDistances, breaks=DISTANCES, plot=FALSE)
empiricalIntraHetDistancesCumulativeFrequencies = cumsum(empiricalIntraHetDistancesHist$counts) / length(empiricalIntraHetDistances)

scatter(list(empiricalIntraHetDistancesHist$mids, DISTANCES), list(empiricalIntraHetDistancesCumulativeFrequencies, freqAtLteDist), "intraHetDistancesDistrib.theoretical_empirical", xlab="Intra-het distance", ylab="Cumulative Frequency", log="x", leg=c("NA12878 HiSeq", "Theoretical"), legPos="topleft", type="b", col=TWO_COLORS)



######################################################################
# Phasing as a function of MEAN intra-het distance:
######################################################################
load("testIntraHetDistances.RData")

hetDistances = list()
phaseRateHetDistances = list()
hetDistancesLeg = vector()

for (nextIndex in 1:length(USE_EMPIRICAL_WINDOWS)) {
	n = USE_EMPIRICAL_WINDOWS[nextIndex]
	meanHetDistNumSitesPhasedConsistentSwitch = scan(paste("~fromer/storage/remove_het_sites.NA12878.HiSeq/meanHetDist_numSites_phased_consistent_switch.chr1.n_", n, ".txt", sep=""), what=list(meanHetDist=0.0, numSites=0, phased=0, consistentPhased=0, switch=0.0))

	hetDistances[nextIndex] = list(meanHetDistNumSitesPhasedConsistentSwitch$meanHetDist)
	phaseRateHetDistances[nextIndex] = list(meanHetDistNumSitesPhasedConsistentSwitch$phased)
	hetDistancesLeg[nextIndex] = paste("Removed hets from HiSeq (window = ", n, ")", sep="")
}

nextIndex = nextIndex+1
hetDistances[nextIndex] = list(MEAN_INTRA_HET_DISTANCES)
phaseRateHetDistances[nextIndex] = list(pPhaseTheta)
hetDistancesLeg[nextIndex] = "Theoretical (window = 2)" # params

scatter(hetDistances, phaseRateHetDistances, "intraHetDistances.theoretical_empirical", xlab="Mean intra-het distance", ylab="Phaseability", leg=hetDistancesLeg, legPos="topright", type="b", col=TWO_COLORS)

scatter(hetDistances, phaseRateHetDistances, "intraHetDistances.log.theoretical_empirical", xlab="Mean intra-het distance", ylab="Phaseability", leg=hetDistancesLeg, legPos="topright", type="b", col=TWO_COLORS, log="y", xlim=c(1, 20000))


######################################################################
# Phasing as a function of window size:
######################################################################
load("theoretical_window_on_empirical.RData")

windows = list()
phaseRateWindows = list()
windowsLeg = vector()

NUM_HET_SITES = as.integer(system("cat ~fromer/storage/phase.NA12878/COMPLETE_LIST.het_sites.interval_list | wc -l", intern=TRUE))
NUM_CHR = as.integer(system("cat ~fromer/storage/phase.NA12878/COMPLETE_LIST.het_sites.interval_list | cut -f1 -d':' | sort | uniq | wc -l", intern=TRUE))
NUM_PHASEABLE_HET_SITES = NUM_HET_SITES - NUM_CHR  # since can't phase the first het site of each chromosome


windowPhasedConsistent = scan(paste("~fromer/storage/phase.NA12878/window_phased_consistent.txt", sep=""), what=list(window=0, phased=0, consistentPhased=0))
windows[1] = list(windowPhasedConsistent$window)
phaseRateWindows[1] = list(windowPhasedConsistent$phased / NUM_PHASEABLE_HET_SITES)
windowsLeg[1] = paste("HiSeq", sep="")


windows[2] = list(WINDOW_SIZES)
phaseRateWindows[2] = list(colMeans(na.omit(phaseProbsPositionWindow)))
windowsLeg[2] = "Theoretical" # params

scatter(windows, phaseRateWindows, "windows.theoretical_empirical", xlab="Window size", ylab="Phaseability", leg=windowsLeg, legPos="topleft", width=14, height=7, type="b", col=TWO_COLORS)



# Use numerical integration over theoretical distances distribution:
load("testWindows.RData")

doneInds = which(pPhaseWindow != -1)

windows[2] = list(WINDOW_SIZES[doneInds])
phaseRateWindows[2] = list(pPhaseWindow[doneInds])
windowsLeg[2] = "Theoretical" # params

scatter(windows, phaseRateWindows, "theoretical_distances.windows.theoretical_empirical", xlab="Window size", ylab="Phaseability", leg=windowsLeg, legPos="topleft", width=14, height=7, type="b", col=TWO_COLORS)



# Use theoretical sampling of distances:
load("theoretical_window.RData")

windows[2] = list(WINDOW_SIZES)
phaseRateWindows[2] = list(colMeans(na.omit(phaseProbsPositionWindow)))
windowsLeg[2] = "Theoretical" # params

scatter(windows, phaseRateWindows, "sampled_theoretical_distances.windows.theoretical_empirical", xlab="Window size", ylab="Phaseability", leg=windowsLeg, legPos="topleft", width=14, height=7, type="b", col=TWO_COLORS)
