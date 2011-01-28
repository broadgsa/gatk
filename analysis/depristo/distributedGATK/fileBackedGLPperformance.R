#d <- read.table("../GATK/trunk/timer.dat", header=T)
require("lattice")
print(xyplot(elapsed.time + delta ~ cycle | name, data=d, scales=list(relation="free"), auto.key=T, type="b", outer=T))