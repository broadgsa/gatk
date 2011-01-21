JOB_START_RATE = 0.1 # chance of starting is 0.1
WORK_UNITS = 100
WORK_RATE = 1
N_TICKS = 300

ticks <- 1:N_TICKS

# the probability that a job starts at exactly tick i
pThreadStartAtTick <- function(i) {
	dexp(i, JOB_START_RATE)
}

jobDoneByI <- function(i) {
	return(sapply(i - ticks, function(x) max(x, 0)) * WORK_RATE)
	#return(pCompleteAtI(i, pStarts, ticks))
}

pThreadDoneByI <- function(i) {
	pStarts <- pThreadStartAtTick(ticks)
	workDoneByThreadStartingAtI <- jobDoneByI(i)
	fracDone <- workDoneByThreadStartingAtI / WORK_UNITS
	doneAtI <- fracDone >= 1
	return(sum(pStarts * doneAtI))
}

pThreadsDoneByI <- function(i, nThreads) {
	pDone <- rep(0, N_TICKS)
	for ( thread : 1:nThreads )
		pDone <- pPrevThreadsNotDoneAtI(pDone, i) + pThreadDoneByI(i)
}

#plot(ticks, workDoneByI(100))
plot(ticks, sapply(ticks, function(i) pThreadDoneByI(i)))

