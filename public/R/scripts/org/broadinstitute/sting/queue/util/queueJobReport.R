library(gsalib)
library(ggplot2)
library(gplots)
library(tools)
library(reshape)

#
# Standard command line switch.  Can we loaded interactively for development
# or executed with RScript
#
args = commandArgs(TRUE)
onCMDLine = ! is.na(args[1])
if ( onCMDLine ) {
  inputFileName = args[1]
  outputPDF = args[2]
} else {
  inputFileName = "Q-26618@gsa4.jobreport.txt"
  #inputFileName = "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/Q-25718@node1149.jobreport.txt"
  #inputFileName = "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/rodPerformanceGoals/history/report.082711.txt"
  outputPDF = NA
}

RUNTIME_UNITS = "(hours)"
ORIGINAL_UNITS_TO_RUNTIME_UNITS = 1/1000/60/60

# 
# Helper function to aggregate all of the jobs in the report across all tables
#
allJobsFromReport <- function(report) {
  names <- c("jobName", "startTime", "analysisName", "doneTime", "exechosts", "runtime")
  sub <- lapply(report, function(table) table[,names])
  do.call("rbind", sub)
}

#
# Creates segmentation plots of time (x) vs. job (y) with segments for the duration of the job
#
plotJobsGantt <- function(gatkReport, sortOverall, includeText) {
  allJobs = allJobsFromReport(gatkReport)
  if ( sortOverall ) {
    title = "All jobs, by analysis, by start time"
    allJobs = allJobs[order(allJobs$analysisName, allJobs$startTime, decreasing=T), ]
  } else {
    title = "All jobs, sorted by start time"
    allJobs = allJobs[order(allJobs$startTime, decreasing=T), ]
  }
  allJobs$index = 1:nrow(allJobs)
  minTime = min(allJobs$startTime)
  allJobs$relStartTime = allJobs$startTime - minTime 
  allJobs$relDoneTime = allJobs$doneTime - minTime
  allJobs$ganttName = paste(allJobs$jobName, "@", allJobs$exechosts)
  maxRelTime = max(allJobs$relDoneTime)
  p <- ggplot(data=allJobs, aes(x=relStartTime, y=index, color=analysisName))
  p <- p + theme_bw()
  p <- p + geom_segment(aes(xend=relDoneTime, yend=index), size=1, arrow=arrow(length = unit(0.1, "cm")))
  if ( includeText )
    p <- p + geom_text(aes(x=relDoneTime, label=ganttName, hjust=-0.2), size=2)
  p <- p + xlim(0, maxRelTime * 1.1)
  p <- p + xlab(paste("Start time, relative to first job", RUNTIME_UNITS))
  p <- p + ylab("Job number")
  p <- p + opts(title=title)
  print(p)
}

#
# Plots scheduling efficiency at job events
#
plotProgressByTime <- function(gatkReport) {
  allJobs = allJobsFromReport(gatkReport)
  nJobs = dim(allJobs)[1]
  allJobs = allJobs[order(allJobs$startTime, decreasing=F),]
  allJobs$index = 1:nrow(allJobs)

  minTime = min(allJobs$startTime)
  allJobs$relStartTime = allJobs$startTime - minTime
  allJobs$relDoneTime = allJobs$doneTime - minTime

  times = sort(c(allJobs$relStartTime, allJobs$relDoneTime))

  countJobs <- function(p) {
    s = allJobs$relStartTime
    e = allJobs$relDoneTime
    x = c() # I wish I knew how to make this work with apply
    for ( time in times )
      x = c(x, sum(p(s, e, time)))
    x
  }

  pending = countJobs(function(s, e, t) s > t)
  done = countJobs(function(s, e, t) e < t)
  running = nJobs - pending - done

  d = data.frame(times=times, pending=pending, running=running, done=done)
  
  p <- ggplot(data=melt(d, id.vars=c("times")), aes(x=times, y=value, color=variable))
  p <- p + facet_grid(variable ~ ., scales="free")
  p <- p + geom_line(size=2)
  p <- p + xlab(paste("Time since start of first job", RUNTIME_UNITS))
  p <- p + opts(title = "Job scheduling")
  print(p)
}

# 
# Creates tables for each job in this group
#
standardColumns = c("jobName", "startTime", "formattedStartTime", "analysisName", "intermediate", "exechosts", "formattedDoneTime", "doneTime", "runtime")
plotGroup <- function(groupTable) {
  name = unique(groupTable$analysisName)[1]
  groupAnnotations = setdiff(names(groupTable), standardColumns)  
  sub = groupTable[,c("jobName", groupAnnotations, "runtime")]
  sub = sub[order(sub$iteration, sub$jobName, decreasing=F), ]
  
  # create a table showing each job and all annotations
  textplot(sub, show.rownames=F)
  title(paste("Job summary for", name, "full itemization"), cex=3)

  # create the table for each combination of values in the group, listing iterations in the columns
  sum = cast(melt(sub, id.vars=groupAnnotations, measure.vars=c("runtime")), ... ~ iteration, fun.aggregate=mean)
  textplot(as.data.frame(sum), show.rownames=F)
  title(paste("Job summary for", name, "itemizing each iteration"), cex=3)

  # histogram of job times by groupAnnotations
  if ( length(groupAnnotations) == 1 && dim(sub)[1] > 1 ) {
    # todo -- how do we group by annotations?
    p <- ggplot(data=sub, aes(x=runtime)) + geom_histogram()
    p <- p + xlab(paste("runtime", RUNTIME_UNITS)) + ylab("No. of jobs")
    p <- p + opts(title=paste("Job runtime histogram for", name))
    print(p)
  }
  
  # as above, but averaging over all iterations
  groupAnnotationsNoIteration = setdiff(groupAnnotations, "iteration")
  if ( dim(sub)[1] > 1 ) {
    try({ # need a try here because we will fail to reduce when there's just a single iteration
      sum = cast(melt(sub, id.vars=groupAnnotationsNoIteration, measure.vars=c("runtime")), ... ~ ., fun.aggregate=c(mean, sd))
      textplot(as.data.frame(sum), show.rownames=F)
      title(paste("Job summary for", name, "averaging over all iterations"), cex=3)
    }, silent=T)
  }
}
    
# print out some useful basic information
print("Report")
print(paste("Project          :", inputFileName))

convertUnits <- function(gatkReportData) {
  convertGroup <- function(g) {
    g$runtime = g$runtime * ORIGINAL_UNITS_TO_RUNTIME_UNITS
    g$startTime = g$startTime * ORIGINAL_UNITS_TO_RUNTIME_UNITS
    g$doneTime = g$doneTime * ORIGINAL_UNITS_TO_RUNTIME_UNITS
    g
  }
  lapply(gatkReportData, convertGroup)
}

#
# Plots runtimes by analysis name and exechosts
#
# Useful to understand the performance of analysis jobs by hosts, 
# and to debug problematic nodes
# 
plotTimeByHost <- function(gatkReportData) {
  fields = c("analysisName", "exechosts", "runtime")
  
  runtimes = data.frame()
  for ( report in gatkReportData ) {
    runtimes = rbind(runtimes, report[,fields])
  }
  
  plotMe <- function(name, vis) {
    p = ggplot(data=runtimes, aes(x=exechosts, y=runtime, group=exechosts, color=exechosts))
    p = p + facet_grid(analysisName ~ ., scale="free")
    p = p + vis()
    p = p + xlab("Job execution host")
    p = p + opts(title = paste(name, "of job runtimes by analysis name and execution host"))
    p = p + ylab(paste("Distribution of runtimes", RUNTIME_UNITS))
    p = p + opts(axis.text.x=theme_text(angle=45, hjust=1, vjust=1))
    print(p)
  }
  
  plotMe("Boxplot", geom_boxplot)
  plotMe("Jittered points", geom_jitter)
}

  
# read the table
gatkReportData <- gsa.read.gatkreport(inputFileName)
gatkReportData <- convertUnits(gatkReportData)
#print(summary(gatkReportData))

if ( ! is.na(outputPDF) ) {
  pdf(outputPDF, height=8.5, width=11)
} 

plotJobsGantt(gatkReportData, T, F)
plotJobsGantt(gatkReportData, F, F)
plotProgressByTime(gatkReportData)
plotTimeByHost(gatkReportData)
for ( group in gatkReportData ) {
  print(group)
 plotGroup(group)
}
  
if ( ! is.na(outputPDF) ) {
  dev.off()
  if (exists("compactPDF")) {
    compactPDF(outputPDF)
  }
}
