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

expected <- function(ps, N, eps) {
	co = 2 * N / ( 1 - eps )
	v = co * ((1 - ps)/(1-eps))^(2 * N - 1)
	v / sum(v)
}

f(afs, c(1,2,3,5,10,50))
x = 1:MAX_AC / 200000
points(x, expected(x,1000,1e-8),type="l",lty=3,lwd=3)

