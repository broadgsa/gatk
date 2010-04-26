#!/usr/bin/env python

import os
import sys
import subprocess
import re
import unittest
import tempfile

def all(iterable):
    for element in iterable:
        if not element:
            return False
    return True

# maximum number of unnamed jobs allowed sa dependencies
MAX_UNNAMED_DEPENDENCIES = 10

class FarmJob:
    def __init__( self, cmd_str_from_user, jobName = None, outputHead = None, outputFile = None, dependencies = [], dependencyNameString = None, dieOnFail = False):
        self.cmd_str_from_user = cmd_str_from_user
        self.jobName = jobName
        self.outputHead = outputHead
        self.outputFile = outputFile
        self.dieOnFail = dieOnFail

        self.dependencies = dependencies
        if self.dependencies == None:
            self.dependencies = []
        elif type(self.dependencies) != list:
            self.dependencies = [self.dependencies]
        self.dependencyNameString = dependencyNameString
        
        if len(self.dependencies) > MAX_UNNAMED_DEPENDENCIES:
            depNames = map(FarmJob.getJobName, self.dependencies)
            if len(filter(None, depNames)) > 1 and len(self.dependencies) != len(filter(None, depNames)):
                # there are some unnamed and some named deps
                raise Exception("Bad job names -- some are named and some are unnamed", depName)
            
        self.jobID = None           # None indicates currently unscheduled
        self.executionString = None # currently unscheduled
        self.executed = False
        self.jobStatus = None
            
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

#    def __repr__(self):
#        return self.__str__()


def longestCommonPrefix(strings):
    #print 'LCP', strings
    if strings == []:
        return ""
    else:
        l = map(len, strings)
        shortestLen = min(l)
        #print '*arg', l, shortestLen
        for i in range(shortestLen):
            c = strings[0][i]
            if not all(map(lambda s: c == s[i], strings)):
                shortestLen = i
                break

        return strings[0][0:shortestLen]

def jobNameWildcard(jobs):
    if len(jobs) == 1:
        return jobs[0].getJobName()
    else:
        return longestCommonPrefix(map(FarmJob.getJobName, jobs)) + "*"

def executeJobs(allJobs, farm_queue = None, just_print_commands = False, debug = True):
    for job in allJobs:
        if type(job) == list:
            # convenience for lists of lists
            map( lambda x: executeJob(x, farm_queue, just_print_commands, debug), job)
        else:
            print 'Preparing to execute', job
            if not job.executed:
                # schedule the dependents
                executeJobs(job.dependencies, farm_queue, just_print_commands, debug = debug)
                executeJob(job, farm_queue, just_print_commands, debug = debug)
                print 'Executed', job

justPrintJobIDCounter = 1

def executeJob(job, farm_queue = None, just_print_commands = False, debug = True, die_on_fail = True):
    global justPrintJobIDCounter
    job.executed = True
    
    # build my execution string
    job.executionString = buildExecutionString(job, farm_queue, debug = debug) 

    if just_print_commands or (globals().has_key("justPrintCommands") and globals().justPrintCommands):
        job.jobID = justPrintJobIDCounter
        justPrintJobIDCounter += 1
    elif farm_queue:
        print 'job.executionString', job.executionString
        result = subprocess.Popen([job.executionString, ""], shell=True, stdout=subprocess.PIPE).communicate()[0]
        p = re.compile('Job <(\d+)> is submitted to queue')
        job.jobID = p.match(result).group(1)
    else:
        # Actually execute the command if we're not just in debugging output mode
        status = os.system(job.executionString)
        if not farm_queue:
            print "<<< Exit code:", status,"\n"
        if die_on_fail != None and status != 0:
            print "### Failed with exit code "+str(status)+" while executing command "+job.cmd_str_from_user
            sys.exit()
        job.jobStatus = int(status)

def buildExecutionString(job, farm_queue = None, debug = True):
    if farm_queue != None:
        if job.outputFile != None:
            farm_stdout = job.outputFile
        elif job.outputHead != None:
            farm_stdout = job.outputHead + ".stdout"
        else:
            #farm_stdout = None
            farm_stdout = "%J.lsf.output"
            
        cmd_str = "bsub -q " + farm_queue 
        if farm_stdout != None:
            cmd_str += " -o " + farm_stdout
        
        # fixme
        if job.dependencies != []:
            cmd_str += buildJobDependencyString(job.dependencies)
        if job.dependencyNameString != None:
            cmd_str += buildJobDependencyString(job.dependencyNameString)
        if job.jobName != None:
            cmd_str += " -J %s" % job.jobName

        cmd_out, cmd_file = tempfile.mkstemp()
        cmd_out = open(cmd_file, 'w')
        cmd_out.write(job.cmd_str_from_user)
        cmd_out.close()

        cmd_str += " < " + cmd_file + ""
        #cmd_str += " '" + job.cmd_str_from_user + "'"
        
        if debug: print ">>> Farming via "+cmd_str
    else:
        cmd_str = job.cmd_str_from_user
        if debug: print ">>> Executing "+cmd_str

    return cmd_str

def allJobsAreNamed(jobs):
    #print map(lambda x: x != None, map(FarmJob.getJobName, jobs))
    return all(map(lambda x: x != None, map(FarmJob.getJobName, jobs)))

def buildJobDependencyString(depJobs):
    if type(depJobs) == str:
        # we are all formally named
        depString = "ended(%s*)" % depJobs
    elif allJobsAreNamed(depJobs):
        # we are all formally named
        depString = "ended(%s)" % jobNameWildcard(depJobs)
    else:
        depString = "&&".join(map(lambda x: "ended(%s)" % x.jobID, depJobs)) 

    return " -w \"%s\"" % depString

def cmd(cmd_str_from_user, farm_queue=False, output_head=None, just_print_commands=False, outputFile = None, waitID = None, jobName = None, die_on_fail = False):
    """if farm_queue != False, submits to queue, other
die_on_fail_msg: if != None, die on command failure (non-zero return) and show die_on_fail_msg"""

    if farm_queue:
        if outputFile <> None:
            farm_stdout = outputFile
        elif output_head <> None:
            farm_stdout = output_head+".stdout"
        else:
            #farm_stdout = None
            farm_stdout = "%J.lsf.output"
            
        cmd_str = "bsub -q "+farm_queue 
        if farm_stdout <> None:
            cmd_str += " -o " + farm_stdout

        if waitID <> None:
            cmd_str += " -w \"ended(%s)\"" % (str(waitID))

        if jobName <> None:
            cmd_str += " -J %s" % (jobName)

        cmd_str += " '"+cmd_str_from_user + "'"
        
        print ">>> Farming via "+cmd_str
    else:
        cmd_str = cmd_str_from_user
        print ">>> Executing "+cmd_str

    if just_print_commands or (globals().has_key("justPrintCommands") and globals().justPrintCommands):
        return -1
    elif farm_queue:
        result = subprocess.Popen([cmd_str, ""], shell=True, stdout=subprocess.PIPE).communicate()[0]
        p = re.compile('Job <(\d+)> is submitted to queue')
        jobid = p.match(result).group(1)
        return jobid        
    else:
        # Actually execute the command if we're not just in debugging output mode
        status = os.system(cmd_str)
        if not farm_queue:
            print "<<< Exit code:", status,"\n"
        if die_on_fail != None and status != 0:
            print "### Failed with exit code "+str(status)+" while executing command "+cmd_str_from_user
            sys.exit()
        return int(status)


# ------------------------------------------------------------------------------------------                                                                                               
#                                                                                                                                                                                          
# Unit testing!                                                                                                                                                                            
#                                                                                                                                                                                          
# ------------------------------------------------------------------------------------------                                                                                               
class TestFarmCommands(unittest.TestCase):
    def setUp(self):
        print ''
        print '-' * 100

    def testMakingJob1(self):
        print 'testMakingJob:'
        job = FarmJob("foo")
        executeJobs([job], just_print_commands = True)

    def testMakingJob2(self):
        print 'testMakingJob2:'
        job = FarmJob("bar", jobName = "barJobs", outputHead = "outputRoot", outputFile = "bar.log", dependencies = [], dieOnFail = True)
        executeJobs([job], "gsa", just_print_commands = True)

    def testDepJobs1(self):
        print 'testDepJobs1:'
        job1 = FarmJob("job1")
        job2 = FarmJob("job2")
        job3 = FarmJob("job3", dependencies = [job1, job2])
        executeJobs([job1, job2, job3], "gsa", just_print_commands = True)

    def testDepJobs2(self):
        print 'testDepJobs2:'
        foojob1 = FarmJob("job1", jobName = "foojob1")
        foojob2 = FarmJob("job2", jobName = "foojob2")
        barjob1 = FarmJob("barjob1", jobName = "barjob1", dependencies = [foojob1])
        barjob2 = FarmJob("barjob2", jobName = "barjob2", dependencies = [foojob1, foojob2])
        bazjob1 = FarmJob("baz1", jobName = "baz1", dependencies = [barjob1, barjob2])
        executeJobs([bazjob1], "gsa", just_print_commands = True)

if __name__ == '__main__':
    unittest.main()

