MAX_AC = 10000
normHist <- function(d, m) {
	x = hist(d$true.ac, breaks=1:20000, plot=F)$counts[1:MAX_AC]
	x / sum(x)
}
	
f <- function(d, acs) {
	cols = rainbow(length(acs), alpha=0.75)
	y = normHist(subset(afs, small.ac == acs[1]))
	x = 1:length(y) / max(d$true.an)
	plot(x, y, type="l", col=cols[1], xlab="True MAF in full population", ylab="Frequency", lwd=3, log="x")
	for (i in 2:length(acs)) {
		points(x, normHist(subset(afs, small.ac == acs[i])), type="l", col=cols[i], lwd=3)
	}
	
	legend("topright", legend=lapply(acs, function(x) paste("AC =", x)), fill=cols, title="Sub-population")
}

expected <- function(maxAN, N, eps, ac1scale = F) {
	scale = 10
	
	f <- function(ps, N) {
		co = 2 * N / ( 1 - eps )
		co * ((1 - ps)/(1-eps))^(2 * N - 1)
	}

	# these are the points that we'll actually show, but we need to do the calculation
	# special for the AC = 1 given the equation actually fits an infinite population
	# not a discrete population with max chromosomes
	ps = 1:maxAN / maxAN
	v = f(ps, N)
	v = v / sum(v)
	
	if ( ac1scale ) {
	  	subps = seq(1, maxAN*scale) / (maxAN * scale)
		#print(subps)
		subv = f(subps, N)
		#print(subv)
		#print(v[1:10])
		pBelowAC1 = sum(subv[1:scale] / sum(subv))
		#print(list(pBelowAC1=pBelowAC1, v1=v[1]))
		v[1] = v[1] + pBelowAC1
	}

	list(ps = ps, pr = v)
}

f(afs, c(1,2,3,5,10,50))

if ( F ) {
scale = 100
ex1 = expected(200000, 1000, 1e-8)
ex2 = expected(200000*scale, 1000, 1e-8)
i = 1:(200000*scale) %% scale == 1
plot(ex2$ps[i], cumsum(ex1$pr), type="l",lty=3,lwd=3, log="x", col="red")
points(ex2$ps[i], cumsum(ex2$pr)[i], type="l",lty=3,lwd=3, log="x")
}

ex = expected(200000, 1000, 1e-8, T)
points(ex$ps, ex$pr, type="l",lty=3,lwd=3)

