args = commandArgs(TRUE);

onCMDLine = ! is.na(args[1])
if (! is.na(args[3]) ) { name = args[3] } else { name = "" }

if ( onCMDLine ) { 
   print(paste("Reading data from", args[1]))
   d = read.table(args[1], header=T, sep="\t")
   d$start.time = as.Date(d$start.time)
   d$end.time = as.Date(d$end.time)
} # only read into d if its' available, otherwise assume the data is already loaded

reportCountingPlot <- function(values, name, moreMargin = 0, ...) {
    par(las=2) # make label text perpendicular to axis
    oldMar <- par("mar")
    par(mar=c(5,8+moreMargin,4,2)) # increase y-axis margin.
    barplot(sort(table(values)), horiz=TRUE, cex.names = 0.5, main = name, xlab="Counts", ...)
    par("mar" = oldMar)
    par("las" = 1)
}

reportHist <- function(values, name, ...) {
    if ( ! all(is.na(values) ) )
        hist(values, main=name, 20, xlab="", col="cornflowerblue", ...)
}

myTable <- function(x, y, reqRowNonZero = F) {
    table <- prop.table(table(x, y), 2)
    ncols = dim(table)[2]

    #print(table)    
    if ( reqRowNonZero )
        table = table[addmargins(table)[1:dim(table)[1],ncols] > 0,]

    return(table)
}

# todo -- must be robust to smaller sizes

plotTable <- function(table, name) {
    ncols = dim(table)[2]
    nrows = dim(table)[1]
    cols = rainbow(nrows)
    tableMin = min(apply(table, 2, min))
    tableMax = max(apply(table, 2, max))
    plot( as.numeric(apply(table, 2, sum)), ylim=c(tableMin, tableMax), type="n", main = name, ylab="Frequency", xlab="Date", xaxt="n")
    axis(1, 1:ncols, labels=colnames(table))
    for ( i in 1:nrows )
        points(table[i,], type="b", col=cols[i])
    legend("topright", row.names(table), fill=cols, cex=0.5)
    #return(table)
}

RUNNING_GATK_RUNTIME <- 60 * 5 #  5 minutes => bad failure

if ( onCMDLine ) pdf(args[2])

successfulRuns <- function(d) {
    x <- rep("Successful", length(d$exception.msg))
    x[! is.na(d$exception.msg)] <- "Failed"
    return(x)
}

generateOneReport <- function(d, header, includeByWeek = T) {
    head <- function(s) {
        return(paste("Section:", header, ":", s))
    }
    
    excepted <- subset(d, exception.msg != "NA")
    badExcepted <- subset(excepted, run.time > RUNNING_GATK_RUNTIME)

    par("mar", c(5, 4, 4, 2))
    frame()
    title(paste("GATK run report", name, "for", Sys.Date(), "\nwith", dim(d)[1], "run repository records"), cex=2)
    
    # cuts by time
    if ( includeByWeek ) {
        plotTable(table(rep("GATK Invocations", length(d$start.time)), cut(d$start.time, "weeks")), head("GATK Invocations by week"))
        plotTable(myTable(successfulRuns(d), cut(d$start.time, "weeks")), head("Successful and failing GATK invocations per week"))
        
        plotTable(myTable(d$svn.version, cut(d$start.time, "weeks")), head("SVN version by week"))
        plotTable(myTable(excepted$walker.name, cut(excepted$start.time, "weeks"), reqRowNonZero = T), head("Walkers with exceptions by week"))
    }
    plotTable(table(rep("GATK Invocations", length(d$start.time)), d$start.time), head("GATK Invocations by day"))
    plotTable(myTable(d$svn.version, d$start.time), head("SVN version by day"))

    reportCountingPlot(d$walker.name, head("Walker invocations"))
    reportCountingPlot(d$svn.version, head("GATK SVN version"))

    # reportCountingPlot(d$java.tmp.directory, head("Java tmp directory"))
    reportCountingPlot(d$working.directory, head("Working directory"))
    reportCountingPlot(d$user.name, head("user"))
    reportCountingPlot(d$host.name, head("host"))
    reportCountingPlot(d$java, head("Java version"))
    reportCountingPlot(d$machine, head("Machine"))
    
    Gb <- 1024^3
    reportHist(d$total.memory / Gb, head("Used memory"))
    reportHist(d$max.memory / Gb, head("Max memory"))
    
    min <- 60
    reportHist(log10(d$run.time / min), head("Run time (log10[min])"))
    
    exceptionColor = "red"
    reportCountingPlot(excepted$walker.name, head("Walker exceptions"), col=exceptionColor)
    reportCountingPlot(subset(excepted, run.time > RUNNING_GATK_RUNTIME)$walker.name, paste(head("Long-running walker exceptions (>"),RUNNING_GATK_RUNTIME,"seconds runtime)"), col=exceptionColor)
    reportCountingPlot(subset(excepted, run.time < RUNNING_GATK_RUNTIME)$walker.name, paste(head("Start-up walker exceptions (<"),RUNNING_GATK_RUNTIME,"seconds runtime)"), col=exceptionColor)
    reportCountingPlot(excepted$user.name, head("Usernames generating exceptions"), col=exceptionColor)
    reportCountingPlot(excepted$exception.msg, head("Exception messages"), 12)
    reportCountingPlot(excepted$exception.at, head("Exception locations"), 12)
}

RUNME = T
if ( RUNME ) {
    lastWeek = levels(cut(d$start.time, "weeks"))[-1]
    generateOneReport(d, "Overall")
    #generateOneReport(subset(d, start.time >= lastWeek), "Just last week to date", includeByWeek = F)
}

if ( onCMDLine ) dev.off()



