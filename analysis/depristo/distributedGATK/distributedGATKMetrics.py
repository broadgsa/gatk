import sys
from optparse import OptionParser
from itertools import *
import random
import re
import datetime

# a simple script that does:
# 1 -- generates a master set of variants following the neutral expectation from a single big population
# 2 -- randomly generates M individuals with variants and genotypes sampled as expected from the big population of variants
# 3 -- writes out the genotypes of these individuals, and their allele frequency
def main():
    global OPTIONS
    usage = "usage: %prog [options] outputFile"
    parser = OptionParser(usage=usage)
    
    (OPTIONS, args) = parser.parse_args()
    if len(args) == 0:
        parser.error("Requires at least one argument")

    print 'file dataset parallel.type nWaysParallel start.time end.time end.to.end.time per.1M.sites job.run.time'
    typere = '.*/(.*).ptype_(\w+).nways_(\d+).*'
    for file in args:
        startTime, endTime, perMSites, runtime = None, None, None, None
        for line in open(file):
            match = re.match(typere, line)
            if match != None: dataset, parallelType, nWays = match.groups()
            startTime = captureStartTime(line, startTime)
            perMSites = capturePerMSites(line, perMSites)
            endTime = captureEndTime(line, endTime)
            runtime = captureRuntime(line, runtime)
        print file, dataset, parallelType, nWays, formatTime(startTime), formatTime(endTime), endToEnd(endTime, startTime), perMSites, runtime
    
def endToEnd(endTime, startTime):
    if endTime < startTime:
        endTime = endTime + datetime.timedelta(1)
    #print 'endToEnd', endTime, startTime
    return total_minutes(endTime - startTime)
    
def formatTime(t):
    return datetime.datetime.strftime(t, formatString)

def total_minutes(td):
    return td.days * 24 * 60 + td.seconds / 60.0

def captureLine(line, regex, func, prevValue):
    match = regex.match(line)
    if match != None:
        if func != None:
            val = func(line)
        else:  
            val = match.group(1)
    else:
        val = None
    #print 'Matching', line, regex, match, prevValue, val

    return val

formatString = "%H:%M:%S"
                
def captureStartTime(line, prev):
    # todo - needs to find the earliest time
    #INFO  11:03:50,202 HelpFormatter - The Genome Analysis Toolkit (GATK) v<unknown>, Compiled <unknown>
    regex = re.compile("INFO\W*(\d+:\d+:\d+).*The Genome Analysis Toolkit.*")
    return selectTime(captureLine(line, regex, None, prev), prev, earlier = True)

def selectTime(newTimeString, oldTime, earlier = False):
    def select():
        if newTimeString == None:
            return oldTime
        else:
            newTime = datetime.datetime.strptime(newTimeString, formatString)
            if oldTime == None:
                return newTime
            elif earlier:
                if newTime < oldTime:
                    return newTime
                else:
                    return oldTime
            else:
                if newTime > oldTime:
                    return newTime
                else:
                    return oldTime
    r = select()
    #if not earlier: print 'selectTime', oldTime, newTimeString, r
    return r


def captureEndTime(line, prev):
    # todo - needs to find the latest time
    regex = re.compile("INFO\W*(\d+:\d+:\d+).*GATKRunReport - Aggregating data for run report.*")
    return selectTime(captureLine(line, regex, None, prev), prev, earlier=False)

unitsToMinutes = {
    'm' : 1.0,
    'h' : 60,
    's' : 1.0/60,
    'd' : 60 * 60
    }

def capturePerMSites(line, prev):
    return captureDoneLine(line, prev, 8, 10)

def captureRuntime(line, prev):
    return captureDoneLine(line, prev, 6, 8)

def captureDoneLine(line, prev, s, e):
        #    INFO  11:04:11,541 TraversalEngine -    chr1:3769010        1.32e+05   20.0 s        2.5 m      1.5%        21.9 m    21.5 m
    regex = re.compile("INFO  .*TraversalEngine -.*done*")
    val = captureLine(line, regex, lambda x: x.split()[s:e], None)
    if val == None: 
        return prev
    else:
        x, u = val
        return float(x) * unitsToMinutes[u]



if __name__ == "__main__":
    main()
