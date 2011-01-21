# todo -- add replicate number to system
# tood -- add scatter gather comparison
d <- read.table("results.dat", header=T)
require("lattice")

subd = data.frame(parallel.type=d$parallel.type, nWaysParallel=d$nWaysParallel, end.to.end.time=d$end.to.end.time,per.1M.sites = d$per.1M.sites, job.run.time = d$job.run.time)

nways = unique(subd$nWaysParallel)
m = max(subd$end.to.end.time)
nNW = subset(subd, end.to.end.time == m)$nWaysParallel[1]
timeAt1 = m * nNW
my.runtime = subset(subd, end.to.end.time == m)$job.run.time[1] * nNW
my.pms = subset(subd, end.to.end.time == m)$per.1M.sites[1]

theo = data.frame(parallel.type="theoretic", end.to.end.time=timeAt1/nways, nWaysParallel=nways, per.1M.sites = my.pms, job.run.time = my.runtime / nways)

subd = rbind(subd, theo)

print(summary(subd))

print(xyplot(end.to.end.time + per.1M.sites + job.run.time ~ nWaysParallel, data=subd[order(subd$nWaysParallel),], group=parallel.type, type="b", outer=T, scale=list(relation="free"), auto.key=T, lwd=c(2,2,1)))

