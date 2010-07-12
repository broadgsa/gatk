import os
import sys
import subprocess
import time
import re
import unittest
import tempfile
import RefseqLibrary

MAX_UNNAMED_DEPENDENCIES = 0
class FarmJob:
    def __init__( self, cmd_str_from_user, jobName, projectName, delayTime = None, jobsDependingOnThis = None,
                  outputFile = None, dependencies = [], dependencyNameString = None, dieOnFail = False, usingFiles = [], memory = None):
        self.cmd_str_from_user = cmd_str_from_user
        self.jobName = jobName
        self.projectName = projectName
        self.jobsDependingOnThis = jobsDependingOnThis
        self.outputFile = outputFile
        self.dieOnFail = dieOnFail
        self.delayTime = delayTime
        self.dependencies = dependencies
        if self.dependencies == None:
            self.dependencies = []
        elif type(self.dependencies) != list:
            self.dependencies = [self.dependencies]
        self.dependencyNameString = dependencyNameString
        self.filesToUse = usingFiles ## provides additional protection for stop-resume spawning by keeping a list of files
        
        if len(self.dependencies) > MAX_UNNAMED_DEPENDENCIES:
            depNames = map(FarmJob.getJobName, self.dependencies)
            if len(filter(None, depNames)) > 1 and len(self.dependencies) != len(filter(None, depNames)):
                # there are some unnamed and some named deps
                raise Exception("Bad job names -- some are named and some are unnamed", depName)
            
        self.jobID = None           # None indicates currently unscheduled
        self.executionString = None # currently unscheduled
        self.executed = False
        self.jobStatus = None
        self.memory = memory.strip("g")

    def __hash__(self):
        return self.cmd_str_from_user.__hash__()
    
    def setDelay(self,delayString):
        self.delayTime = delayString

    def getJobName(self): 
        return self.jobName

    def getJobIDString(self): 
        if self.jobName == None:
            if self.jobID == None:
                return "UNNAMED"
            else:
                return str(self.jobID)
        else:
            return self.getJobName()
    
    def __str__(self):
        return "[JOB: name=%s id=%s depending on (%s) with cmd=%s]" % (self.getJobName(), self.jobID, ','.join(map(FarmJob.getJobIDString, self.dependencies)), self.cmd_str_from_user)

    def getNumberOfDependencies(self):
        return len(self.dependencies)

    def getNumberOfDependers(self):
        if ( self.jobsDependingOnThis == None ):
            return 0
        else:
            return len(jobsDependingOnThis)

def compareDependence(job1,job2):
    if ( job1.getJobName() in job2.dependencies ):
        return 1
    elif ( job2.getJobName() in job1.dependencies ):
        return -1
    else:
        return job1.getNumberOfDependers() - job2.getNumberOfDependers()

class JobDispatchError(Exception):

    def __init__(self,value):
        self.value = value

    def __str__(self):
        return repr(self.value)

def delayToGlobalTime(delayStr):
    # delayStr is of the form 1:2:3 for 1 day 2 hrs 3 minutes
    dayHrMin = delayStr.split(":")
    additionalSecs = 24*60*60*int(dayHrMin[0]) + 60*60*int(dayHrMin[1]) + 60*int(dayHrMin[2])
    futureTime = time.localtime(time.time()+additionalSecs)
    return ":".join([str(futureTime[0]),str(futureTime[1]),str(futureTime[2]),str(futureTime[3]),str(futureTime[4])])

def buildDependencyString(jobList):
    dep = '\''
    paddedList = list()
    for name in jobList:
        paddedList.append('ended("'+name+'")')
    dep += " && ".join(paddedList)+'\''
    return dep

def buildSubmitString(farmJob,queue):
    submitStr = "bsub -q "+queue
    if ( farmJob.jobName != None ):
        submitStr += " -J "+farmJob.jobName
    if ( farmJob.projectName != None ):
        submitStr += " -P "+farmJob.projectName
    if ( farmJob.outputFile != None ):
        submitStr += " -o "+farmJob.outputFile
    if ( farmJob.delayTime != None ):
        submitStr += " -b "+delayToGlobalTime(farmJob.delayTime)
    if ( farmJob.dependencies != [] ):
        submitStr += " -w "+buildDependencyString(farmJob.dependencies)
    if ( farmJob.memory != None ):
        submitStr += " -R \"rusage[mem="+farmJob.memory+"]\""
    submitStr += " "+farmJob.cmd_str_from_user
    return submitStr

def writeResumeFile(jobList,start,max,filepath,specialHash = None):
    print("Writing start was: "+str(start+max))
    try:
        output = open(filepath,'w')
        printToFile = True
    except IOError:
        print("Unable to open resume file, "+filepath+" dumping to stdout")
        printToFile = False
    hashVal = 0
    if ( specialHash != None ):
        hashFunction = specialHash
    else:
        hashFunction = hash
    for i in range(start,len(jobList)):
        hashVal = hashVal ^ hashFunction(jobList[i])
    if ( printToFile ):
        output.write(str(start+max)+"\n")
        output.write(str(hashVal)+"\n")
        output.write(jobList[start].cmd_str_from_user)
    else:
        print(str(start))
        print(str(hashVal))
        print(jobList[start].cmd_str_from_user)

def writeResumeFileFinal(filepath):
    output = open(filepath,'w')
    output.write("ALL_JOBS_HAVE_BEEN_SPAWNED")
    output.close()

def checkResumeFile(filepath,jobs,specialHash = None):
    if ( filepath == None or not os.path.exists(filepath) ):
        return 0
    else:
        input = open(filepath)
        firstLine = input.readline()
        if ( firstLine == "ALL_JOBS_HAVE_BEEN_SPAWNED" ):
            raise JobDispatchError("All jobs for this project were spawned.")
        else:
            start = int(firstLine)
        print("reading start was "+str(start))
        hashValToMatch = int(input.readline())
        input.close()
        if ( specialHash == None ):
            hashFunction = hash
        else:
            hashFunction = specialHash
        hashVal = 0
        for i in range(start,len(jobs)):
            hashVal = hashVal ^ hashFunction(jobs[i])
        if ( hashVal == hashValToMatch ):
            return start
        else:
            str1="The commands for remaining jobs hashed to "+str(hashVal)+" but previous hash was "+str(hashValToMatch)+"."
            str2="Please check the resume file "+filepath+" to see that the job command has not changed."
            raise JobDispatchError(str1+" "+str2)

def hashJobAndIntervals(farmJob):
    fjhash = farmJob.__hash__()
    interval_file = farmJob.filesToUse
    intervals = list()
    for line in open(interval_file[0]):
        if ( not line.startswith("@") ):
            spline = line.strip().split()
            try:
                intervals.append(Interval(spline[0],int(spline[1]),int(spline[2])))
            except IndexError:
                print(line)
                raise IndexError("List index out of range")
    inthash = 0
    for interval in intervals:
        inthash = inthash ^ interval.__hash__()
    return fjhash ^ inthash
    
class JobDispatcher:
    def __init__(self,lsf_queues = ["long"], queue_limits = dict([["long",500]]), exceed_total_limit_action = "fail", print_only = False, action_string = None):
        self.queues = lsf_queues
        self.limits = queue_limits
        self.action = exceed_total_limit_action
        self.action_string = action_string
        self.spawned = list()
        if ( len(queue_limits) > 0 ):
            self.maxJobs = sum(queue_limits.values())
        else:
            self.maxJobs = -1
        self.print_only = print_only
        self.startDays = 0
        self.startHours = 0
        self.startMins = 0

    def setInitialDelay(self,delay):
        dhm = delay.split(":")
        self.startDays = int(dhm[0])
        self.startHours = int(dhm[1])
        self.startMins = int(dhm[2])

    # external accessor, sorts by dependency and ensures user hasn't exceeded the set limits
    def dispatchAll(self,farmJobs):
        farmJobs.sort(compareDependence)
        if ( self.maxJobs > -1 and len(farmJobs) > self.maxJobs ):
            if ( self.action == "space" and self.action_string != None):
                self._dispatchWithSpacing(farmJobs)
            elif ( self.action == "resume" and self.action_string != None ):
                self._dispatchWithStopResume(farmJobs)
            else:
                raise JobDispatchError("Number of jobs to dispatch, "+str(len(farmJobs))+", exceeds maximum ("+str(self.maxJobs)+").")
        else:
            self._dispatchAll(farmJobs)

    # internal accessor, loops over queues and dispatches jobs up to the limit
    def _dispatchAll(self,farmJobs):
        for queue in self.queues:
            dispatchedToQueue = 0
            while ( dispatchedToQueue < self.limits[queue] and len(farmJobs) > 0 ):
                farmJob = farmJobs.pop(0)
                self._dispatch(farmJob,queue)
                self.spawned.append(farmJob)
                dispatchedToQueue = 1 + dispatchedToQueue

    # internal dispatch accessor; job limits dealt with at this point
    def _dispatch(self,job,queue):
        if ( self.print_only ):
            print(buildSubmitString(job,queue))
        else:
            lsf_response = subprocess.Popen([buildSubmitString(job,queue), ""], shell=True, stdout=subprocess.PIPE).communicate()[0]
            p = re.compile('Job <(\d+)> is submitted to queue')
            job.jobID = p.match(lsf_response).group(1)
            print(lsf_response)

    # spaces jobs out by using the delay command to LSF
    def _dispatchWithSpacing(self,farmJobs):
        dayHrMin = self.action_string.split(":")
        days = self.startDays
        hours = self.startHours
        mins = self.startMins
        dispatchBuffer = list()
        while ( len(farmJobs) > 0 ):
            nextJob = farmJobs.pop(0)
            nextJob.setDelay(":".join([str(days),str(hours),str(mins)]))
            dispatchBuffer.append(nextJob)
            if ( len(dispatchBuffer) >= self.maxJobs or len(farmJobs) == 0 ):
                self._dispatchAll(dispatchBuffer)
                days += int(dayHrMin[0])
                hours += int(dayHrMin[1])
                mins += int(dayHrMin[2])
                dispatchBuffer = list()

    # dispatches jobs up to the limit, outputs a file describing where to resume again
    # if the file already exists, resumes at that job
    # file also contains a hashcode for the remaining jobs -- if these don't match the
    # hashcode calculated from farmJobs, dispatching is aborted via an exception
    def _dispatchWithStopResume(self,farmJobs,jobHash = None):
        startAt = checkResumeFile(self.action_string,farmJobs,jobHash)
        dispatchBuffer = list()
        while ( len(dispatchBuffer) < self.maxJobs and len(farmJobs) > startAt ):
            dispatchBuffer.append(farmJobs.pop(startAt))
        self._dispatchAll(dispatchBuffer)
        if ( len(farmJobs) > startAt ):
            writeResumeFile(farmJobs,startAt,self.maxJobs,self.action_string,jobHash)
        else:
            writeResumeFileFinal(self.action_string)
class Interval:
    def __init__(self,chrom,start,stop):
        self.chromosome = chrom
        self.start = start
        self.stop = stop

    def size(self):
        return self.stop - self.start

    def bedFormat(self):
        return "\t".join([self.chromosome,str(self.start),str(self.stop),"+","target_whatever"])

    def __str__(self):
        return self.chromosome+":"+str(self.start)+"-"+str(self.stop)

    def __hash__(self):
        return self.__str__().__hash__()

class GATKDispatcher(JobDispatcher):

    def __init__(self,jarfile,memory,walker,args,output_directory,reference = None, bams = None, intervals = None,
                 queues = ["long"], limits = dict([["long",500]]), print_only = False, action = "fail", delay = None):
        self.jarfile = jarfile
        self.memory = memory
        self.walker = walker
        self.args = args
        self.reference = reference
        self.bams = bams
        self.intervals = intervals
        self.outputDir = output_directory
        self.project = "GSA_GATK_Analysis"
        if ( action == "resume" ):
            action_string = self.outputDir+"GATKDispatcher/resumeJobs.txt"
        elif ( action == "space" ):
            action_string = delay
        else:
            action_string = None
        JobDispatcher.__init__(self,queues,limits,action,print_only,action_string)
        self.check()
        self.baseCommand = "java -Xmx"+self.memory+" -jar "+self.jarfile+" -T "+self.walker+" "+args

    def check(self):
        if ( not os.path.exists(self.jarfile) ):
            raise JobDispatchError("The provided GATK jarfile "+str(self.jarfile)+" does not exist")
        if ( self.intervals != None and not os.path.exists(self.intervals) ):
            raise JobDispatchError("The provided interval list file, "+str(self.intervals)+" does not exist.")
        if ( self.bams != None and not os.path.exists(self.bams) ):
            raise JobDispatchError("The provided bam, or bam list, "+str(self.bams)+" does not exist.")
        if ( self.reference == None ):
            print("Warning: No reference supplied to GATKDispatcher. Most analyses require a reference.")

    def setProject(self,newProject):
        self.project = newProject

    def addReadFilter(self,filter,filter_args = None):
        if ( filter_args == None ):
            self.baseCommand += " -rf "+filter
        else:
            raise JobDispatchError("GATK Dispatcher does not yet support read filters with arguments (e.g. blacklisting)")

    # remove bam files, reference, intervals, read filters, etc
    def resetCommand(self):
        self.baseCommand = "java -Xmx"+self.memory+" -jar "+self.jarfile+" -T "+self.walker+" "+args

    def dispatchByInterval(self,base_limit):
        self.check()
        dispatchCommand = self.baseCommand + " -R "+self.reference
        if ( self.bams != None ):
            dispatchCommand += " -I "+self.bams
        intervals = open(self.intervals)
        job_number = 0
        bases_for_job = 0
        intervals_for_job = list()
        headerLines = list()
        farmJobs = list()
        if ( not os.path.exists(self.outputDir+"GATKDispatcher/") ):
            os.mkdir(self.outputDir+"GATKDispatcher/")
        for line in intervals:
            if ( line.startswith("@") ):
                headerLines.append(line)
            else:
                spline = line.strip().split()
                chrom = spline[0]
                start = int(spline[1])
                stop = int(spline[2])
                interval = Interval(chrom,start,stop)
                intervals_for_job.append(interval)
                bases_for_job += interval.size()
                if ( bases_for_job > base_limit ):
                    farmJobs.append(self._buildIntervalJob(job_number,headerLines,intervals_for_job,dispatchCommand))
                    intervals_for_job = list()
                    bases_for_job = 0
                    job_number += 1
        if ( len(intervals_for_job) > 0 ): ## there's still some leftover
            farmJobs.append(self._buildIntervalJob(job_number,headerLines,intervals_for_job,dispatchCommand))
        self.dispatchAll_Interval(farmJobs)

    def dispatchAll_Interval(self,jobs):
        if ( self.action == "resume" ):
            self._dispatchWithStopResume(jobs,hashJobAndIntervals)
        else:
            self.dispatchAll(jobs)

    def _buildIntervalJob(self,num,header,intervals,command):
        job_dir = self.outputDir+"GATKDispatcher/dispatch"+str(num)+"/"
        cmd = self.appendOutput(command,num,job_dir)
        if ( not os.path.exists(job_dir) ):
            os.mkdir(job_dir)
        intFile = open(job_dir+"job"+str(num)+"_intervals.interval_list",'w')
        intFile.write("".join(header))
        for interval in intervals:
            intFile.write(interval.bedFormat()+"\n")
        intFile.close()
        cmd += " -L "+job_dir+"job"+str(num)+"_intervals.interval_list"
        job = FarmJob(cmd,self.project+"_job"+str(num),self.project,None,None,job_dir+"bsub_out.txt",list(),None,False,[job_dir+"job"+str(num)+"_intervals.interval_list"],self.memory)
        return job

    def appendOutput(self,command,num,job_dir):
        if ( self.walker == "UnifiedGenotyper" ):
            cmdToDispatch = command + " -varout "+job_dir+"job"+str(num)+"_calls.vcf"
        elif ( self.walker == "CoverageStatistics" or self.walker == "DepthOfCoverage" ):
            cmdToDispatch = command + " -o "+job_dir+"job"+str(num)
        elif ( self.walker == "CombineDuplicates" or self.walker == "TableRecalibration" or self.walker == "ClipReads" ):
            cmdToDispatch = command + " -o "+job_dir+"job"+str(num)+"_output.bam"
        else:
            cmdToDispatch = command + " -o "+job_dir+"job"+str(num)+".txt"
        return cmdToDispatch

    def dispatchByGene(self, geneNames):
        self.genes = RefseqLibrary.getRefseqGenes(geneNames)
        dispatchCommand = self.baseCommand + " -R "+self.reference
        farmJobs = list()
        jobNumber = ""
        headerLines = RefseqLibrary.getIntervalHeaderLines()

        if ( not os.path.exists(self.outputDir+"GATKDispatcher/") ):
            os.mkdir(self.outputDir+"GATKDispatcher/")
        if ( self.bams != None ):
            dispatchCommand += " -I "+self.bams

        for gene in self.genes:
            jobNumber = "_"+gene.getGeneName()
            intervals = gene.getExonIntervals()
            farmJobs.append(self._buildIntervalJob(jobNumber,headerLines,intervals,dispatchCommand))
            
        self.dispatchAll_Interval(farmJobs)

    def dispatchByTargetDesign(self,designFile):
        self.genes = RefseqLibrary.parseDesignFile(designFile)
        dispatchCommand = self.baseCommand + " -R "+self.reference
        farmJobs = list()
        jobNumber = ""
        headerLines = RefseqLibrary.getIntervalHeaderLines()

        if ( not os.path.exists(self.outputDir+"GATKDispatcher/") ):
            os.mkdir(self.outputDir+"GATKDispatcher/")
        if ( self.bams != None ):
            dispatchCommand += " -I "+self.bams

        for gene in self.genes:
            jobNumber = "_"+gene.getGeneName()
            intervals = gene.getExonIntervals()
            farmJobs.append(self._buildIntervalJob(jobNumber,headerLines,intervals,dispatchCommand))

        self.dispatchAll_Interval(farmJobs)
