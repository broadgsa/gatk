import os.path
import sys
from optparse import OptionParser
from itertools import *
from xml.etree.ElementTree import *
import gzip

MISSING_VALUE = "NA"
RUN_REPORT_LIST = "GATK-run-reports"

def main():
    global OPTIONS
    usage = "usage: %prog [options] mode file1 ... fileN"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", dest="verbose",
                        action='store_true', default=False,
                        help="If provided, verbose progress will be enabled")         

    parser.add_option("-o", "--o", dest="output",
                        type='string', default=None,
                        help="if provided, output will go here instead of stdout")
         
    (OPTIONS, args) = parser.parse_args()
    if len(args) == 0:
        parser.error("Requires at least GATKRunReport xml to analyze")

    stage = args[0]
    files = resolveFiles(args[1:])

    # open up the output file
    if OPTIONS.output != None:
        out = openFile(OPTIONS.output,'w')
    else:
        out = sys.stdout

    handler = getHandler(stage)(stage, out)
    handler.initialize(files)

    # parse all of the incoming files
    for report in readReports(files):
        # todo -- add matching here
        handler.processRecord(report)

    handler.finalize(files)
    out.close()

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
    
# def 
class RecordAsTable(StageHandler):
    def __init__(self, name, out):
        StageHandler.__init__(self, name, out)
        
    def initialize(self, args):
        self.fields = list()
        self.formatters = dict()
    
        def id(elt): return elt.text
        def toString(elt): return '"%s"' % elt.text
        
        def formatExceptionMsg(elt):
            return '"%s"' % elt.find("message").text

        def formatExceptionAt(elt):
            return '"%s"' % elt.find("stacktrace").find("string").text
        
        def add(names, func):
            for name in names:
                addComplex(name, [name], [func])

        def addComplex(key, fields, funcs):
            self.fields.extend(fields)
            self.formatters[key] = zip(fields, funcs)
    
        add(["id", "walker-name", "svn-version"], id)
        add(["start-time", "end-time"], toString)      
        add(["run-time", "java-tmp-directory", "working-directory", "user-name", "host-name"], id)
        add(["java", "machine"], toString)
        add(["max-memory", "total-memory", "iterations", "reads"], id)
        addComplex("exception", ["exception-msg", "exception-at"], [formatExceptionMsg, formatExceptionAt])
        # add(["command-line"], toString)    
        
        print >> self.out, "\t".join(self.fields)

    def processRecord(self, record):
        parsed = parseReport(record, self.formatters)

        def oneField(field):
            val = MISSING_VALUE
            if field in parsed:
                val = parsed[field]
            return val

        print >> self.out, "\t".join([ oneField(field) for field in self.fields ])

def parseReport(report, allFormatters):
    bindings = dict()
    for elt in report:
        if elt.tag in allFormatters:
            fieldFormats = allFormatters[elt.tag]
            # we actually care about this tag
            for field, formatter in fieldFormats:
                bindings[field] = formatter(elt)
    return bindings


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
        for arg in args:
            print 'Deleting file: ', arg
            os.remove(arg)

addHandler('archive', Archive)


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

def readReports(files):
    #print files
    for file in files:
        input = openFile(file)
        tree = ElementTree(file=input)
        elem = tree.getroot()
        if elem.tag == RUN_REPORT_LIST:
            for sub in elem:
                yield sub
        else:
            yield elem
    
if __name__ == "__main__":
    main()