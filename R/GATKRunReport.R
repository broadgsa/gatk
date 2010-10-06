args = commandArgs(TRUE);

onCMDLine = ! is.na(args[1])
if (! is.na(args[3]) ) { name = args[3] } else { name = "" }

if ( onCMDLine ) { 
   print(paste("Reading data from", args[1]))
   d = read.table(args[1], header=T, sep="\t")
   #d$start.time = as.Date(d$start.time)
   d$end.time = as.Date(d$end.time)
} # only read into d if its' available, otherwise assume the data is already loaded

reportCountingPlot <- function(values, name, moreMargin = 0, ...) {
    par(las=2) # make label text perpendicular to axis
    oldMar <- par("mar")
    par(mar=c(5,8+moreMargin,4,2)) # increase y-axis margin.
    barplot(sort(table(factor(values))), horiz=TRUE, cex.names = 0.5, main = name, xlab="Counts", log="x", ...)
    par("mar" = oldMar)
    par("las" = 1)
}

reportConditionalCountingPlot <- function(values, conditions, name, moreMargin = 0, ...) {
    par(las=2) # make label text perpendicular to axis
    oldMar <- par("mar")
    par(mar=c(5,8+moreMargin,4,2)) # increase y-axis margin.
    t = table(values, conditions)
    t = t[, order(colSums(t))]
    #print(list(t = t))
    nconds = dim(t)[2]
    cols = rainbow(nconds)
    barplot(t, legend.text = T, horiz=TRUE, cex.names = 0.5, main = name, xlab="Counts", col=cols, cex=0.5, ...)
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

plotTable <- function(table, name, ...) {
    ncols = dim(table)[2]
    nrows = dim(table)[1]
    if ( ! is.null(nrows) ) {
        cols = rainbow(nrows)
        tableMin = min(apply(table, 2, min))
        tableMax = max(apply(table, 2, max))
        plot( as.numeric(apply(table, 2, sum)), ylim=c(tableMin, tableMax), type="n", main = name, ylab="Frequency", xlab="Date", xaxt="n", ...)
        axis(1, 1:ncols, labels=colnames(table))
        for ( i in 1:nrows )
            points(table[i,], type="b", col=cols[i])
        legend("topright", row.names(table), fill=cols, cex=0.5)
        #return(table)
    }
}

RUNNING_GATK_RUNTIME <- 60 * 5 #  5 minutes => bad failure

if ( onCMDLine ) pdf(args[2])

successfulRuns <- function(d) {
    x <- rep("Successful", length(d$exception.msg))
    x[d$exception.msg != "NA" & d$is.user.exception == "true"] <- "Failed with UserException"
    x[d$exception.msg != "NA" & d$is.user.exception == "false"] <- "Failed with StingException"
    x[d$exception.msg != "NA" & (d$is.user.exception == "NA" | is.na(d$is.user.exception))] <- "Failed with StingException before UserException code"
    return(x)
}

addSection <- function(name) {
    par("mar", c(5, 4, 4, 2))
    frame()
    title(name, cex=2)
}

dropit <- function (d, columns = names(d), ...)
{
   d[columns] = lapply(d[columns], "[", drop=TRUE, ...)
   d
}

generateOneReport <- function(d, header, includeByWeek = T) {
    head <- function(s) {
        return(paste("Section:", header, "\n", s))
    }
    
    excepted <- dropit(subset(d, exception.msg != "NA"))
    UserExceptions <- dropit(subset(excepted, is.user.exception == "true"))
    StingExceptions <- dropit(subset(excepted, is.user.exception == "false" | is.user.exception == "NA" | is.na(is.user.exception)))

    addSection(paste("GATK run report", name, "for", Sys.Date(), "\nwith", dim(d)[1], "run repository records"))

    reportCountingPlot(d$walker.name, head("Walker invocations"))
    reportConditionalCountingPlot(d$user.name, d$walker.name, head("Walker invocations by user"))
    reportCountingPlot(d$svn.version, head("SVN version"))
	reportConditionalCountingPlot(d$svn.version, d$user.name, head("SVN by user"))

    
    # cuts by time
    if ( includeByWeek ) {
        plotTable(table(rep("GATK Invocations", length(d$end.time)), cut(d$end.time, "weeks")), head("GATK Invocations by week"))
        plotTable(myTable(successfulRuns(d), cut(d$end.time, "weeks")), head("Successful and failing GATK invocations per week"))
        
        plotTable(myTable(d$svn.version, cut(d$end.time, "weeks")), head("SVN version by week"))
    }
    plotTable(table(rep("GATK Invocations", length(d$end.time)), d$end.time), head("GATK Invocations by day"))
    plotTable(myTable(d$svn.version, d$end.time), head("SVN version by day"))

	# 
	# Exception handling
	#
	addExceptionSection <- function(subd, subname, exceptionColor) { 	
		addSection(paste(subname))
		#print(list(subd = length(subd$end.time), name=subname))
		reportCountingPlot(subd$walker.name, head(paste("Walkers with", subname)), col=exceptionColor)
    	reportCountingPlot(subd$exception.at, head(paste(subname, "locations")), 12, col=exceptionColor)
    	#reportCountingPlot(subd$exception.msg, head(paste(subname, "messages")), 12, col=exceptionColor)
		reportConditionalCountingPlot(subd$user.name, subd$exception.at, head(paste("Walker invocations by user for", subname)), 12)

	    if ( includeByWeek && length(subd$end.time) > 0 ) {
    	    plotTable(myTable(subd$walker.name, cut(subd$end.time, "weeks"), reqRowNonZero = T), head(paste("Walkers with", subname,"by week")), col=exceptionColor)
		}
	}

	addExceptionSection(excepted, "Exceptions", "grey")
    reportCountingPlot(excepted$user.name, head("Usernames generating exceptions"), col="grey")

	addExceptionSection(StingExceptions, "StingExceptions", "red")
	addExceptionSection(UserExceptions, "UserExceptions", "blue")

    
    Gb <- 1024^3
    reportHist(d$total.memory / Gb, head("Used memory"))
    reportHist(d$max.memory / Gb, head("Max memory"))
    
    min <- 60
    reportHist(log10(d$run.time / min), head("Run time (log10[min])"))

    reportCountingPlot(d$user.name, head("user"))
    reportCountingPlot(d$host.name, head("host"))

    reportCountingPlot(d$java, head("Java version"))
    reportCountingPlot(d$machine, head("Machine"))
    reportCountingPlot(d$working.directory, head("Working directory"))
}

RUNME = T
if ( RUNME ) {
    lastWeek = levels(cut(d$end.time, "weeks"))[-1]
    generateOneReport(d, "Overall")
    #generateOneReport(subset(d, end.time >= lastWeek), "Just last week to date", includeByWeek = F)
}

if ( onCMDLine ) dev.off()



