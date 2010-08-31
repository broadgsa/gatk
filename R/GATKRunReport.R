args = commandArgs(TRUE);

onCMDLine = ! is.na(args[1])

if ( onCMDLine ) { 
   print(paste("Reading data from", args[1]))
   d = read.table(args[1], header=T, sep="\t")
} # only read into d if its' available, otherwise assume the data is already loaded

reportCountingPlot <- function(values, name, moreMargin = 0, ...) {
    par(las=2) # make label text perpendicular to axis
    par(mar=c(5,8+moreMargin,4,2)) # increase y-axis margin.
    barplot(sort(table(values)), horiz=TRUE, cex.names = 0.5, main = name, xlab="Counts", ...)
}

reportHist <- function(values, name, ...) {
    hist(values, main=name, 20, xlab="", col="cornflowerblue", ...)
}

RUNNING_GATK_RUNTIME <- 60 * 5 #  5 minutes => bad failure
excepted <- subset(d, exception.msg != "NA")
badExcepted <- subset(excepted, run.time > RUNNING_GATK_RUNTIME)

if ( onCMDLine ) pdf(args[2])
reportCountingPlot(d$walker.name, "Walker invocations")
reportCountingPlot(d$svn.version, "GATK SVN version")
reportCountingPlot(d$java.tmp.directory, "Java tmp directory")
reportCountingPlot(d$working.directory, "Working directory")
reportCountingPlot(d$user.name, "User")
reportCountingPlot(d$host.name, "host")
reportCountingPlot(d$java, "Java version")
reportCountingPlot(d$machine, "Machine")

Gb <- 1024^3
reportHist(d$total.memory / Gb, "Used memory")
reportHist(d$max.memory / Gb, "Max memory")

min <- 60
reportHist(log10(d$run.time / min), "Run time (log10[min])")

exceptionColor = "red"
reportCountingPlot(excepted$walker.name, "Walker exceptions", col=exceptionColor)
reportCountingPlot(subset(excepted, run.time > RUNNING_GATK_RUNTIME)$walker.name, paste("Long-running walker exceptions (>",RUNNING_GATK_RUNTIME,"seconds runtime)"), col=exceptionColor)
reportCountingPlot(subset(excepted, run.time < RUNNING_GATK_RUNTIME)$walker.name, paste("Start-up walker exceptions (<",RUNNING_GATK_RUNTIME,"seconds runtime)"), col=exceptionColor)
reportCountingPlot(excepted$user.name, "Usernames generating exceptions", col=exceptionColor)
reportCountingPlot(excepted$exception.msg, "Exception messages", 12)
reportCountingPlot(excepted$exception.at, "Exception locations", 12)

if ( onCMDLine ) dev.off()



