from __future__ import with_statement
import farm_commands
import os.path
import sys
from optparse import OptionParser
import picard_utils
from gatkConfigParser import *
import re
from itertools import *
import math

def phredQScore( nMismatches, nBases ):
    #print 'phredQScore', nMismatches, nBases
    if nMismatches == 0:
        return 40
    elif nBases == 0:
        return 0
    else:
        return -10 * math.log10(float(nMismatches) / nBases)
        

def phredScore2ErrorProp(qual):
    #print 'phredScore2ErrorProp', qual
    return math.pow(10.0, float(qual) / -10.0)

expectedHeader = 'rg,dn,Qrep,pos,NBases,MMismatches,Qemp'.split(',')
defaultValues = '0,**,0,0,0,0,0'.split(',')
class RecalData(dict):

    def __init__(self):
        self.parse(expectedHeader, defaultValues)

    def parse(self, header, data):
        # rg,dn,Qrep,pos,NBases,MMismatches,Qemp
        types = [str, str, int, int, int, int, int]
        for head, expected, datum, type in zip(header, expectedHeader, data, types):
            if head <> expected:
                raise ("Unexpected header in rawData %s %s %s" % (head, expected, datum))
            #print 'Binding => ', head, type(datum)
            self[head] = type(datum)
        #print self
        return self
    
    def set(self, header, values):
        for head, val in zip(header, values):
            self[head] = val
            
    def __getattr__(self, name):
        return self[name]
    
    # rg,dn,Qrep,pos,NBases,MMismatches,Qemp
    def readGroup(self): return self.rg
    def dinuc(self): return self.dn
    def qReported(self): return self.Qrep
    def qEmpirical(self): return self.Qemp
    def cycle(self): return self.pos
    def nBases(self): return self.NBases
    def nMismatches(self): return self.MMismatches
    def nExpectedMismatches(self): return self.nBases() * phredScore2ErrorProp(self.qReported())

    def combine(self, moreData):
        # grab useful info
        sumErrors = self.nExpectedMismatches()
        for datum in moreData:
            self.NBases += datum.nBases()
            self.MMismatches += datum.nMismatches()
            sumErrors += datum.nExpectedMismatches()
        self.updateQemp()
        self.Qrep = phredQScore(sumErrors, self.nBases())
        #print 'self.Qrep is now', self.Qrep
        return self
        
    def updateQemp(self):
        newQemp = phredQScore( self.nMismatches(), self.nBases() )
        #print 'Updating qEmp', self.Qemp, newQemp
        self.Qemp = newQemp
        return newQemp

    def __str__(self):
        return "rg=%s cycle=%d dinuc=%s qrep=%.1f qemp=%.1f nbases=%d nmismatchs=%d" % ( self.readGroup(), self.cycle(), self.dinuc(), self.qReported(), self.qEmpirical(), self.nBases(), self.nMismatches())
    
#    def __init__(dinuc, Qrep, pos, nbases, nmismatches, qemp ):
#        self.dinuc = dinuc
#        self.Qrep = Qrep

def rawDataStream(file):
    """Yields successive lists containing the CSVs in the data file; excludes headers"""
    header = None
    for line in open(file):
        if line.find("#") <> -1: continue
        else:
            data = line.strip().split(',')
            if line.find("rg,") <> -1:
                header = data
            else:
                yield RecalData().parse(header, data)

def rawDataByReadGroup(rawDataFile):
    """Yields a stream of the data in rawDataFile, grouped by readGroup"""
    for readGroup, generator in groupby(rawDataStream(rawDataFile), key=RecalData.readGroup):
        yield (readGroup, list(generator))

def combineRecalData(separateData):
    return RecalData().combine(separateData)

def groupRecalData(allData, key=None):
    s = sorted(allData, key=key)
    values = [ [key, combineRecalData(vals)] for key, vals in groupby(s, key=key) ]
    return sorted( values, key=lambda x: x[0])

#
# let's actually analyze the data!
#
def analyzeReadGroup(readGroup, data, outputRoot):
    print 'Read group => ', readGroup
    print 'Number of elements => ', len(data)

    files = []
    if OPTIONS.toStdout:
        qReportedVsqEmpirical(readGroup, data, sys.stdout )
        qDiffByCycle(readGroup, data, sys.stdout)
        qDiffByDinuc(readGroup, data, sys.stdout)
    else:
        def outputFile(tail):
            file = outputRoot + tail
            files.append(file)
            return file
            
        with open(outputFile(".empirical_v_reported_quality.dat"), 'w') as output:
            qReportedVsqEmpirical(readGroup, data, output )
        with open(outputFile(".quality_difference_v_cycle.dat"), 'w') as output:
            qDiffByCycle(readGroup, data, output)
        with open(outputFile(".quality_difference_v_dinucleotide.dat"), 'w') as output:
            qDiffByDinuc(readGroup, data, output)
            
    print 'Files', files
    return analyzeFiles(files)

def qDiffByCycle(readGroup, allData, output):
    print >> output, '# Note Qreported is a float here due to combining Qreported across quality bins -- Qreported is the expected Q across all Q bins, weighted by nBases'
    print >> output, 'Cycle    Qreported   Qempirical     Qempirical_Qreported     nMismatches     nBases'
    for cycle, datum in groupRecalData(allData, key=RecalData.cycle):
        datum.set(['rg', 'dn', 'pos'], [readGroup, '**', cycle])
        diff = datum.qEmpirical() - datum.qReported()
        print >> output, "%d   %2.2f  %2.2f   %2.2f     %12d    %12d" % (datum.cycle(), datum.qReported(), datum.qEmpirical(), diff, datum.nMismatches(), datum.nBases())

# 
#     public void qualityDiffVsDinucleotide(PrintStream file, final String readGroup) {
#         ArrayList<RecalData> ByCycle = new ArrayList<RecalData>();
#         ArrayList<MeanReportedQuality> ByCycleReportedQ = new ArrayList<MeanReportedQuality>();
#         file.printf("dinuc,Qemp-obs,Qemp,Qobs,B,N%n");
#         RecalData All = new RecalData(0,0,readGroup,"");
#         MeanReportedQuality AllReported = new MeanReportedQuality();
#         for (int c=0; c < RecalData.NDINUCS; c++) {
#             ByCycle.add(new RecalData(-1, -1,readGroup,RecalData.dinucIndex2bases(c)));
#             ByCycleReportedQ.add(new MeanReportedQuality());
#         }
# 
#         for ( RecalData datum: getRecalData(readGroup) ) {
#             int dinucIndex = RecalData.dinucIndex(datum.dinuc); //bases2dinucIndex(datum.dinuc.charAt(0), datum.dinuc.charAt(1), false);
#             ByCycle.get(dinucIndex).inc(datum.N, datum.B);
#             ByCycleReportedQ.get(dinucIndex).inc(datum.qual, datum.N);
#             All.inc(datum.N, datum.B);
#             AllReported.inc(datum.qual, datum.N);
#         }
# 
#         for (int c=0; c < RecalData.NDINUCS; c++) {
#             double empiricalQual = -10 * Math.log10((double)ByCycle.get(c).B / ByCycle.get(c).N);
#             double reportedQual = ByCycleReportedQ.get(c).result();
#             file.printf("%s, %f, %f, %f, %d, %d%n", ByCycle.get(c).dinuc, empiricalQual-reportedQual, empiricalQual, reportedQual, ByCycle.get(c).B, ByCycle.get(c).N);
#         }
#     }

def qDiffByDinuc(readGroup, allData, output):
    print >> output, '# Note Qreported is a float here due to combining Qreported across quality bins -- Qreported is the expected Q across all Q bins, weighted by nBases'
    print >> output, 'Dinuc    Qreported   Qempirical     Qempirical_Qreported     nMismatches     nBases'
    for dinuc, datum in groupRecalData(allData, key=RecalData.dinuc):
        datum.set(['rg', 'dn', 'pos'], [readGroup, dinuc, '*'])
        diff = datum.qEmpirical() - datum.qReported()
        print >> output, "%s   %2.2f  %2.2f   %2.2f     %12d    %12d" % (datum.dinuc(), datum.qReported(), datum.qEmpirical(), diff, datum.nMismatches(), datum.nBases())

def qReportedVsqEmpirical(readGroup, allData, output):
    print >> output, 'Qreported    Qempirical   nMismatches     nBases'
    for key, datum in groupRecalData(allData, key=RecalData.qReported):
        datum.set(['rg', 'dn', 'Qrep', 'pos'], [readGroup, '**', key, '*'])
        print >> output, "%2d  %2.2f   %12d    %12d" % (datum.qReported(), datum.qEmpirical(), datum.nMismatches(), datum.nBases())

def analyzeRawData(rawDataFile):
    for readGroup, data in rawDataByReadGroup(rawDataFile):
        if OPTIONS.selectedReadGroups == [] or readGroup in OPTIONS.selectedReadGroups:
            root, sourceFilename = os.path.split(rawDataFile)
            if ( OPTIONS.outputDir ): root = OPTIONS.outputDir
            outputRoot = os.path.join(root, "%s.%s.%s" % ( sourceFilename, readGroup, 'analysis' ))
            analyzeReadGroup(readGroup, data, outputRoot)

plottersByFile = {
    "raw_data.csv$" : analyzeRawData,
    "recal_data.csv$" : analyzeRawData,
    "empirical_v_reported_quality" : 'PlotQEmpStated',
    "quality_difference_v_dinucleotide" : 'PlotQDiffByDinuc',
    "quality_difference_v_cycle" : 'PlotQDiffByCycle' }

def getPlotterForFile(file):
    for pat, analysis in plottersByFile.iteritems():
        if re.search(pat, file):
            if type(analysis) == str:
                return config.getOption('R', analysis, 'input_file')
            else:
                analysis(file)
    return None

def analyzeFiles(files):
    #print 'analyzeFiles', files
    Rscript = config.getOption('R', 'Rscript', 'input_file')
    for file in files:
        plotter = getPlotterForFile(file)
        if plotter <> None:
            cmd = ' '.join([Rscript, plotter, file])
            farm_commands.cmd(cmd, None, None, just_print_commands = OPTIONS.dry)

if __name__ == "__main__":
    usage = """usage: %prog -c config.cfg files*"""

    parser = OptionParser(usage=usage)
    parser.add_option("-q", "--farm", dest="farmQueue",
                        type="string", default=None,
                        help="Farm queue to send processing jobs to")
    parser.add_option("-d", "--dir", dest="outputDir",
                        type="string", default=None,
                        help="If provided, analysis output files will be written to this directory")
    parser.add_option("-c", "--config", dest="configs",
                        action="append", type="string", default=[],
                        help="Configuration file")                        
    parser.add_option("-s", "--stdout", dest="toStdout",
                        action='store_true', default=False,
                        help="If provided, writes output to standard output, not to files")
    parser.add_option("", "--dry", dest="dry",
                        action='store_true', default=False,
                        help="If provided, nothing actually gets run, just a dry run")
    parser.add_option("-g", "--readGroup", dest="selectedReadGroups",
                        action="append", type="string", default=[],
                        help="If provided, only the provided read groups will be analyzed")                        
                        
    (OPTIONS, args) = parser.parse_args()
    #if len(args) != 3:
    #    parser.error("incorrect number of arguments")

    if len(OPTIONS.configs) == 0:
        parser.error("Requires at least one configuration file be provided")

    config = gatkConfigParser(OPTIONS.configs)
  
    analyzeFiles(args)