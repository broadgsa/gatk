plot1 <- function(d, name) {
	d = subset(d, dataset == name)
	subd = data.frame(parallel.type=d$parallel.type, nWaysParallel=d$nWaysParallel, end.to.end.time=d$end.to.end.time,per.1M.sites = d$per.1M.sites, job.run.time = d$job.run.time)

	nways = unique(subd$nWaysParallel)
	m = max(subset(subd, nWaysParallel == min(nways))$end.to.end.time)
	nNW = subset(subd, end.to.end.time == m)$nWaysParallel[1]
	timeAt1 = m * nNW
	my.runtime = subset(subd, end.to.end.time == m)$job.run.time[1] * nNW
	my.pms = subset(subd, end.to.end.time == m)$per.1M.sites[1]

	theo = data.frame(parallel.type="theoretic", end.to.end.time=timeAt1/nways, nWaysParallel=nways, per.1M.sites = my.pms, job.run.time = my.runtime / nways)

	subd = rbind(subd, theo)

	print(summary(subd))

	print(xyplot(log10(end.to.end.time) + per.1M.sites + log10(job.run.time) ~ log2(nWaysParallel), data=subd[order(subd$nWaysParallel),], group=parallel.type, type="b", outer=T, scale=list(relation="free"), auto.key=T, lwd=c(2,2,1), main=name))
	
	return(subd)
}

myData <- read.table("results.new.dat", header=T)
require("lattice")

for (name in unique(d$dataset))
	plot1(myData, name)
