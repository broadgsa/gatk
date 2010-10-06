import os.path
import sys
from optparse import OptionParser
from itertools import *
from xml.etree.ElementTree import *
import gzip
import datetime
import re

MISSING_VALUE = "NA"
RUN_REPORT_LIST = "GATK-run-reports"

def main():
    global OPTIONS
    usage = "usage: %prog [options] mode file1 ... fileN"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v", "--verbose", dest="verbose",
                        action='store_true', default=False,
                        help="If provided, verbose progress will be enabled")         

    parser.add_option("", "--overwrite", dest="overwrite",
                        action='store_true', default=False,
                        help="If provided, archive mode will overwrite destination file, if it exists (DANGEROUS)")         

    parser.add_option("-o", "--o", dest="output",
                        type='string', default=None,
                        help="if provided, output will go here instead of stdout")

    parser.add_option("", "--no-dev", dest="noDev",
                        action='store_true', default=False,
                        help="if provided, only records not coming from a dev version of GATK will be included")

    parser.add_option("", "--stingExceptionsOnly", dest="stingExceptionsOnly",
                        action='store_true', default=False,
                        help="if provided, only records with StingExceptions will be generated in report")

    parser.add_option("", "--max_days", dest="maxDays",
                        type='int', default=None,
                        help="if provided, only records generated within X days of today will be included")

    parser.add_option("-D", "--delete_while_archiving", dest="reallyDeleteInArchiveMode",
                        action='store_true', default=False,
                        help="if provided, we'll actually delete records when running in archive mode")
         
    (OPTIONS, args) = parser.parse_args()
    if len(args) == 0:
        parser.error("Requires at least GATKRunReport xml to analyze")

    stage = args[0]
    files = resolveFiles(args[1:])

    # open up the output file
    if OPTIONS.output != None:
        if stage == "archive" and os.path.exists(OPTIONS.output) and not OPTIONS.overwrite:
            raise "archive output file already exists, aborting!", OPTIONS.output
        out = openFile(OPTIONS.output,'w')
    else:
        out = sys.stdout

    handler = getHandler(stage)(stage, out)
    handler.initialize(files)

    # parse all of the incoming files
    counter = 0
    for report in readReports(files):
        # todo -- add matching here
        handler.processRecord(report)
        counter += 1

    handler.finalize(files)
    if OPTIONS.output != None: out.close()
    print 'Processed records:', counter 

#
# Stage HANDLERS
#
class StageHandler:
    def __init__(self, name, out):
        self.name = name
        self.out = out
        
    def getName(self): return self.name
        
    def initialize(self, args):
        pass # print 'initialize'
        
    def processRecord(self, record):
        pass # print 'processing record', record
        
    def finalize(self, args):
        pass # print 'Finalize'


# a map from stage strings -> function to handle record
HANDLERS = dict()
def addHandler(name, handler):
    HANDLERS[name] = handler
    
def getHandler(stage):
    return HANDLERS[stage]

def eltIsException(elt):
    return elt.tag == "exception"
    
def parseException(elt):
    msgElt = elt.find("message")
    msgText = "MISSING"
    userException = "NA"
    if msgElt != None: msgText = msgElt.text
    stackTrace = elt.find("stacktrace").find("string").text
    if elt.find("is-user-exception") != None:
    	#print elt.find("is-user-exception")
        userException = elt.find("is-user-exception").text
    return msgText, stackTrace, userException

def javaExceptionFile(javaException):
	m = re.search("\((.*\.java:.*)\)", javaException)
	if m != None:
		return m.group(1)
	else:
		javaException

class RecordDecoder:
    def __init__(self):
        self.fields = list()
        self.formatters = dict()
    
        def id(elt): return elt.text
        def toString(elt): return '%s' % elt.text
        
        def formatExceptionMsg(elt):
            return '%s' % parseException(elt)[0]

        def formatExceptionAt(elt):
            return '%s' % parseException(elt)[1]

        def formatExceptionAtBrief(elt):
            return '%s' % javaExceptionFile(parseException(elt)[1])

        def formatExceptionUser(elt):
            return '%s' % parseException(elt)[2]
        
        def add(names, func):
            for name in names:
                addComplex(name, [name], [func])

        def addComplex(key, fields, funcs):
            self.fields.extend(fields)
            self.formatters[key] = zip(fields, funcs)
    
        add(["id", "walker-name", "svn-version", "phone-home-type"], id)
        add(["start-time", "end-time"], toString)      
        add(["run-time", "java-tmp-directory", "working-directory", "user-name", "host-name"], id)
        add(["java", "machine"], toString)
        add(["max-memory", "total-memory", "iterations", "reads"], id)
        addComplex("exception", ["exception-msg", "exception-at", "exception-at-brief", "is-user-exception"], [formatExceptionMsg, formatExceptionAt, formatExceptionAtBrief, formatExceptionUser])
        # add(["command-line"], toString)          
    
    def decode(self, report):
        bindings = dict()
        for elt in report:
            if elt.tag in self.formatters:
                fieldFormats = self.formatters[elt.tag]
                # we actually care about this tag
                for field, formatter in fieldFormats:
                    bindings[field] = formatter(elt)

        # add missing data
        for field in self.fields:
            if field not in bindings:
                bindings[field] = MISSING_VALUE

        return bindings
    
# def 
class RecordAsTable(StageHandler):
    def __init__(self, name, out):
        StageHandler.__init__(self, name, out)
        
    def initialize(self, args):
        self.decoder = RecordDecoder()
        print >> self.out, "\t".join(self.decoder.fields)

    def processRecord(self, record):
        parsed = self.decoder.decode(record)

        def oneField(field):
            val = MISSING_VALUE
            if field in parsed:
                val = parsed[field]
                if val.find(" ") != -1:
                    val = "\"" + val + "\""
            return val

        print >> self.out, "\t".join([ oneField(field) for field in self.decoder.fields ])

addHandler('table', RecordAsTable)

class RecordAsXML(StageHandler):
    def __init__(self, name, out):
        StageHandler.__init__(self, name, out)

    def initialize(self, args):
        print >> self.out, "<%s>" % RUN_REPORT_LIST
        
    def processRecord(self, record):
        print >> self.out, tostring(record)

    def finalize(self, args):
        print >> self.out, "</%s>" % RUN_REPORT_LIST

addHandler('xml', RecordAsXML)

class Archive(RecordAsXML):
    def __init__(self, name, out):
        RecordAsXML.__init__(self, name, out)

    def finalize(self, args):
        RecordAsXML.finalize(self, args)
        for arg in args:
            if OPTIONS.verbose: print 'Deleting file: ', arg
            if OPTIONS.reallyDeleteInArchiveMode:
                os.remove(arg)
        print 'Deleted', len(args), 'files'
        
addHandler('archive', Archive)

class ExceptionReport(StageHandler):
    #FIELDS = ["Msg", "At", "SVN.versions", "Walkers", 'Occurrences', 'IDs']
    def __init__(self, name, out):
        StageHandler.__init__(self, name, out)
        self.exceptions = []

    def initialize(self, args):
        self.decoder = RecordDecoder()
        #print >> self.out, "\t".join(self.FIELDS)
        
    def processRecord(self, record):
        for elt in record:
            if eltIsException(elt):
                self.exceptions.append(self.decoder.decode(record))
                break

    def finalize(self, args):
        commonExceptions = list()
        
        def addToCommons(ex):
            for common in commonExceptions:
                if common.equals(ex):
                    common.update(ex)
                    return
            commonExceptions.append(CommonException(ex))

        for ex in self.exceptions:
            addToCommons(ex)
        commonExceptions = sorted(commonExceptions, None, lambda x: x.counts)   
            
        for common in commonExceptions:
            msg, at, svns, walkers, counts, ids, duration, users, userError = common.toStrings()
            
            if OPTIONS.stingExceptionsOnly and userError == "true":
                continue
                
            print >> self.out, ''.join(['*'] * 80)
            print >> self.out, 'Exception              :', msg
            print >> self.out, '    is-user-exception? :', userError
            print >> self.out, '    at                 :', at
            print >> self.out, '    walkers            :', walkers
            print >> self.out, '    svns               :', svns
            print >> self.out, '    duration           :', duration
            print >> self.out, '    occurrences        :', counts
            print >> self.out, '    users              :', users
            print >> self.out, '    ids                :', ids
            
class CommonException:
    MAX_SET_ITEMS_TO_SHOW = 5

    def __init__(self, ex):
        self.msgs = set([ex['exception-msg']])
        self.at = ex['exception-at']
        self.svns = set([ex['svn-version']])
        self.users = set([ex['user-name']])
        self.userError = ex['is-user-exception']
        self.counts = 1
        self.times = set([decodeTime(ex['end-time'])])
        self.walkers = set([ex['walker-name']])
        self.ids = set([ex['id']])
        
    def equals(self, ex):
        return self.at == ex['exception-at']
        
    def update(self, ex):
        self.msgs.add(ex['exception-msg'])
        self.svns.add(ex['svn-version'])
        self.users.add(ex['user-name'])
        self.counts += 1
        self.walkers.add(ex['walker-name'])
        self.times.add(decodeTime(ex['end-time']))
        self.ids.add(ex['id'])

    def bestExample(self, examples):
        def takeShorter(x, y):
            if len(y) < len(x):
                return y
            else:
                return x
        return reduce(takeShorter, examples)
        
    def setString(self, s):
        if len(s) > self.MAX_SET_ITEMS_TO_SHOW:
            s = [x for x in s][0:self.MAX_SET_ITEMS_TO_SHOW] + ["..."]
        return ','.join(s)
        
    def duration(self):
        x = sorted(filter(lambda x: x != "ND", self.times))
        if len(x) >= 2:
            return "-".join(map(lambda x: x.strftime("%m/%d/%y"), [x[0], x[-1]]))
        elif len(x) == 1:
            return x[0]
        else:
            return "ND"
            
        
    def toStrings(self):
        return [self.bestExample(self.msgs), self.at, self.setString(self.svns), self.setString(self.walkers), self.counts, self.setString(self.ids), self.duration(), self.setString(self.users), self.userError] 

addHandler('exceptions', ExceptionReport)

  
  
class SummaryReport(StageHandler):
    #FIELDS = ["Msg", "At", "SVN.versions", "Walkers", 'Occurrences', 'IDs']
    def __init__(self, name, out):
        StageHandler.__init__(self, name, out)
        self.reports = []

    def initialize(self, args):
        self.decoder = RecordDecoder()
        #print >> self.out, "\t".join(self.FIELDS)
        
    def processRecord(self, record):
        self.reports.append(self.decoder.decode(record))

    def finalize(self, args):
        print >> self.out, 'GATK run summary for          :', datetime.datetime.today()
        print >> self.out, '    number of runs            :', len(self.reports)
        print >> self.out, '    number of StingExceptions :', len(filter(isStingException, self.reports))
        print >> self.out, '    number of UserExceptions  :', len(filter(isUserException, self.reports))
        print >> self.out, '    users                     :', ', '.join(set(map(userID, self.reports)))

def userID(rec):
    return rec['user-name']

def isStingException(rec):
    return rec['exception-at'] != "NA" and rec['is-user-exception'] == "false"

def isUserException(rec):
    return rec['exception-at'] != "NA" and rec['is-user-exception'] == "true"

addHandler('summary', SummaryReport)  
  
#
# utilities
#
def openFile(filename, mode='r'):
    if ( filename.endswith(".gz") ):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def resolveFiles(paths):
    allFiles = list()
    def resolve1(path):
        if not os.path.exists(path):
            raise Exception("Path doesn't exist: " + path)
        elif os.path.isfile(path):
            allFiles.append(path)
        else:
            def one(arg, dirname, files):
                #print dirname, files
                #print dirname
                allFiles.extend(map( lambda x: os.path.join(path, x), files ))
                #print files
    
            os.path.walk(path, one, None)

    map( resolve1, paths )
    return allFiles

def decodeTime(time):
    if time == "ND":
        return "ND"
    else:
        return datetime.datetime.strptime(time.split()[0], "%Y/%m/%d")
    #return datetime.datetime.strptime(time, "%Y/%m/%d %H.%M.%S")

def passesFilters(elt):
    if OPTIONS.noDev and eltTagEquals(elt,'build-type','dev'):
        return False
    if OPTIONS.maxDays != None:
        now = datetime.datetime.today()
        now = datetime.datetime(now.year, now.month, now.day)
        #    <start-time>2010/08/31 15.38.00</start-time>
        eltTime = decodeTime(elt.find('end-time').text)
        diff = now - eltTime
        #print eltTime, now, diff, diff.days
        if diff.days > OPTIONS.maxDays:
            return False

    return True
    
def readReports(files):
    #print files
    for file in files:
        if OPTIONS.verbose: print 'Reading file', file
        input = openFile(file)
        try:
            tree = ElementTree(file=input)
        except:
            print "Ignoring excepting file", file
            continue

        elem = tree.getroot()
        if elem.tag == RUN_REPORT_LIST:
            for sub in elem:
                if passesFilters(sub):
                    yield sub
        else:
            if passesFilters(elem):
                yield elem
    
if __name__ == "__main__":
    main()
