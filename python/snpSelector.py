import os.path
import sys
from optparse import OptionParser
from vcfReader import *
#import pylab
from itertools import *
import math
import random
    
DEBUG = False

class CallCovariate:
    def __init__(self, feature, left, right, FPRate = None, cumulative = False):
        self.feature = feature
        self.left = left

        if cumulative:
            self.right = '*'
        else:
            self.right = right

        self.FPRate = FPRate
        
    def containsVariant(self, call):
        fieldVal = call.getField(self.feature)
        return fieldVal >= self.left and (self.right == '*' or fieldVal <= self.right)

    def getFPRate(self): return self.FPRate
    def getFeature(self): return self.feature
    
    def getCovariateField(self): return self.getFeature() + '_RQ'
    
    def __str__(self): return "[CC feature=%s left=%s right=%s]" % (self.feature, self.left, self.right)

class RecalibratedCall:
    def __init__(self, call, features):
        self.call = call
        self.features = dict([[feature, None] for feature in features])
        
    def recalFeature( self, feature, FPRate ):
        assert self.features[feature] == None, "Feature " + feature + ' has value ' + str(self.features[feature]) + ' for call ' + str(self.call) # not reassigning values
        assert FPRate <= 1 and FPRate >= 0
        self.features[feature] = FPRate
        
    def getFeature( self, feature, missingValue = None, phredScaleValue = False ):
        v = self.features[feature]
        if v == None:
            return missingValue
        elif phredScaleValue:
            return phredScale(v)
        else:
            return v
        
    def jointFPErrorRate(self):
        #print self.features
        logTPRates = [math.log10(1-r) for r in self.features.itervalues() if r <> None]
        logJointTPRate = reduce(lambda x, y: x + y, logTPRates, 0)
        logJointTPRate = min(logJointTPRate, 1e-3 / 3) # approximation from het of 0.001
        jointTPRate = math.pow(10, logJointTPRate)
        #print logTPRates
        #print logJointTPRate, jointTPRate
        return 1 - jointTPRate
        
    def featureStringList(self):
        return ','.join(map(lambda feature: '%s=Q%d' % (feature, self.getFeature(feature, '*', True)), self.features.iterkeys()))        
        
    def __str__(self):
        return '[%s: %s => Q%d]' % (str(self.call), self.featureStringList(), phredScale(self.jointFPErrorRate()))

def readVariants( file, maxRecords = None, decodeAll = True, downsampleFraction = 1 ):
    f = open(file)
    header, ignore, lines = readVCFHeader(f)

    def parseVariant(args):
        header1, VCF, counter = args
        if random.random() <= downsampleFraction:
            return VCF
        else:
            return None

    return header, ifilter(None, imap(parseVariant, islice(lines2VCF(lines, extendedOutput = True, decodeAll = decodeAll), maxRecords)))

def selectVariants( variants, selector = None ):
    if selector <> None:
        return filter(selector, variants)
    else:
        return variants

def titv(variants):
    ti = len(filter(VCFRecord.isTransition, variants))
    tv = len(variants) - ti
    titv = ti / (1.0*max(tv,1))

    return titv

def dbSNPRate(variants):
    inDBSNP = len(filter(VCFRecord.isKnown, variants))
    return float(inDBSNP) / len(variants)

def gaussian(x, mu, sigma):    
    constant = 1 / math.sqrt(2 * math.pi * sigma**2)
    exponent = -1 * ( x - mu )**2 / (2 * sigma**2)
    return constant * math.exp(exponent)

# if target = T, and FP calls have ti/tv = 0.5, we want to know how many FP calls
# there are in N calls with ti/tv of X.  
# 
def titvFPRateEstimate(variants, target):
    titvRatio = titv(variants)
    
    # f <- function(To,T) { (To - T) / (1/2 - T) + 0.001 }
    def theoreticalCalc():
        if titvRatio >= target:
            FPRate = 0
        else:
            FPRate = (titvRatio - target) / (0.5 - target)
        FPRate = min(max(FPRate, 0), 1)
        TPRate = max(min(1 - FPRate, 1 - dephredScale(OPTIONS.maxQScore)), dephredScale(OPTIONS.maxQScore))
        if DEBUG: print 'FPRate', FPRate, titvRatio, target
        assert FPRate >= 0 and FPRate <= 1
        return TPRate
    
    # gaussian model
    def gaussianModel():
        LEFT_HANDED = True
        sigma = 5
        constant = 1 / math.sqrt(2 * math.pi * sigma**2)
        exponent = -1 * ( titvRatio - target )**2 / (2 * sigma**2)
        TPRate = gaussian(titvRatio, target, sigma) / gaussian(target, target, sigma)
        if LEFT_HANDED and titvRatio >= target:
            TPRate = 1
        TPRate -= dephredScale(OPTIONS.maxQScore)
        if DEBUG: print 'TPRate', TPRate, constant, exponent, dephredScale(OPTIONS.maxQScore)
        return TPRate
    
    #denom = (0.2 - 0.8 * titvRatio)
    #FPRate = 1
    #if denom <> 0:
    #    FPRate = (1.0 / (target+1)) * (titvRatio - target) / denom

    FPRate = 1 - gaussianModel()
    nVariants = len(variants)
    
    if DEBUG: print ':::', nVariants, titvRatio, target, FPRate
    
    return titvRatio, FPRate
    
def phredScale(errorRate):
    return -10 * math.log10(max(errorRate, 1e-10))

def dephredScale(qscore):
    return math.pow(10, float(qscore) / -10)

def frange6(*args):
    """A float range generator."""
    start = 0.0
    step = 1.0

    l = len(args)
    if l == 1:
        end = args[0]
    elif l == 2:
        start, end = args
    elif l == 3:
        start, end, step = args
        if step == 0.0:
            raise ValueError, "step must not be zero"
    else:
        raise TypeError, "frange expects 1-3 arguments, got %d" % l

    v = start
    while True:
        if (step > 0 and v >= end) or (step < 0 and v <= end):
            raise StopIteration
        yield v
        v += step

def compareFieldValues( v1, v2 ):
    if type(v1) <> type(v2):
        #print 'Different types', type(v1), type(v2)
        c = cmp(type(v1), type(v2))
    else:
        c = cmp(v1, v2)
    #print 'Comparing %s %s = %s' % (v1, v2, c)
    return c

def calculateBins(variants, field, minValue, maxValue, partitions):
    sortedVariants = sorted(variants, key = lambda x: x.getField(field))  # cmp = compareFieldValues, 
    sortedValues = map(lambda x: x.getField(field), sortedVariants)
    
    targetBinSize = len(variants) / (1.0*partitions)
    #print sortedValues
    uniqBins = groupby(sortedValues)
    binsAndSizes = map(lambda x: [x[0], len(list(x[1]))], uniqBins)
    #print 'BINS AND SIZES', binsAndSizes

    def bin2Break(bin): return [bin[0], bin[0], bin[1]]
    bins = [bin2Break(binsAndSizes[0])]
    for bin in binsAndSizes[1:]:
        #print '  Breaks', bins
        #print '  current bin', bin
        curSize = bin[1]
        prevSize = bins[-1][2]
        #print curSize, prevSize
        if curSize + prevSize > targetBinSize:
            #print '     => appending', bin2Break(bin)
            bins.append(bin2Break(bin))
        else:
            bins[-1][1] = bin[0]
            bins[-1][2] += curSize

    #print 'Returning ', bins
    #sys.exit(1)
    return bins

def fieldRange(variants, field):
    values = map(lambda v: v.getField(field), variants)
    minValue = min(values)
    maxValue = max(values)
    #rangeValue = maxValue - minValue
    bins = calculateBins(variants, field, minValue, maxValue, OPTIONS.partitions)
    validateBins(bins)
    return minValue, maxValue, bins

def validateBins(bins):
    #print 'Bins are', bins
    for left1, right1, count1 in bins:
        for left2, right2, count2 in bins:
            def contains2(x):
                return left2 < x and x < right2

            if left1 <> left2 and right1 <> right2:
                if None in [left1, left2, right1, right2]:
                    pass # we're ok
                elif contains2(left1) or contains2(right2):
                    raise Exception("Bad bins", left1, right1, left2, right2)

def printFieldQualHeader(more = ""):
    print '  field left right nVariants nNovels titv titvNovels dbSNP fprate q', more
    
def printFieldQual( field, left, right, variants, FPRate, more = ""):
    novels = selectVariants(variants, VCFRecord.isNovel)
    print '  %s %s %8d %8d %.2f %.2f %.2f %.2e %d' % (field, binString(left, right), len(variants), len(novels), titv(variants), titv(novels), dbSNPRate(variants), FPRate, phredScale(FPRate)), more

def binString(left, right):
    leftStr = str(left)
    if type(left) == float: leftStr = "%.2f" % left
    rightStr = "%5s" % str(right)
    if type(right) == float: rightStr = "%.2f" % right
    return '%8s %8s' % (leftStr, rightStr)


#
#
#
def recalibrateCalls(variants, fields, callCovariates):
    def phred(v): return int(round(phredScale(v)))
    
    for variant in variants:
        recalCall = RecalibratedCall(variant, fields) 
        originalQual = variant.getField('QUAL') 

        for callCovariate in callCovariates:
            if callCovariate.containsVariant(variant):
                FPR = callCovariate.getFPRate()
                recalCall.recalFeature(callCovariate.getFeature(), FPR)
                recalCall.call.setField(callCovariate.getCovariateField(), phred(FPR))

                
        recalCall.call.setField('QUAL', phred(recalCall.jointFPErrorRate())) 
        recalCall.call.setField('OQ', originalQual)
        #print 'recalibrating', variant.getLoc()
        #print '  =>',  variant
        yield recalCall.call
    
#
#
#
def optimizeCalls(variants, fields, titvTarget):
    callCovariates = calibrateFeatures(variants, fields, titvTarget)
    recalCalls = recalibrateCalls(variants, fields, callCovariates)
    return recalCalls, callCovariates

def printCallQuals(recalCalls, titvTarget, info = ""):
    print '--------------------------------------------------------------------------------'
    print info
    calibrateFeatures(recalCalls, ['QUAL'], titvTarget, printCall = True, cumulative = False, forcePrint = True )
    print 'Cumulative'
    calibrateFeatures(recalCalls, ['QUAL'], titvTarget, printCall = True, cumulative = True, forcePrint = True )



def all( p, l ):
    for elt in l:
        if not p(elt): return False
    return True

def variantBinsForField(variants, field):
    #if not all( lambda x: x.hasField(field), variants):
    #    raise Exception('Unknown field ' + field)
    
    minValue, maxValue, bins = fieldRange(variants, field)
    if DEBUG: print 'Field range', minValue, maxValue
    if DEBUG: print 'Partitions', bins
    return bins

def mapVariantBins(variants, field, cumulative = False):
    bins = variantBinsForField(variants, field)
    
    def variantsInBin(bin):
        cc = CallCovariate(field, bin[0], bin[1], cumulative = cumulative)

        return cc.left, cc.right, selectVariants(variants, lambda v: cc.containsVariant(v))
        
    return imap( variantsInBin, bins )

def calibrateFeatures(variants, fields, titvTarget, printCall = False, cumulative = False, forcePrint = False ):
    covariates = []    

    printFieldQualHeader()
    for field in fields:
        if DEBUG: print 'Optimizing field', field
        
        titv, FPRate = titvFPRateEstimate(variants, titvTarget)
        #print 'Overall FRRate:', FPRate, nErrors, phredScale(FPRate)

        for left, right, selectedVariants in mapVariantBins(variants, field, cumulative = cumulative):
            if len(selectedVariants) > max(OPTIONS.minVariantsPerBin,1) or forcePrint:
                titv, FPRate = titvFPRateEstimate(selectedVariants, titvTarget)
                dbsnp = dbSNPRate(selectedVariants)
                covariates.append(CallCovariate(field, left, right, FPRate))
                printFieldQual(field, left, right, selectedVariants, FPRate )
            else:
                print 'Not calibrating bin', left, right, 'because it contains too few variants:', len(selectedVariants)

    return covariates

class CallCmp:
    def __init__(self, nTP, nFP, nFN):
        self.nTP = nTP
        self.nFP = nFP
        self.nFN = nFN
    
    def FPRate(self):
        return (1.0*self.nFP) / max(self.nTP + self.nFP, 1)

    def FNRate(self):
        return (1.0*self.nFN) / max(self.nTP + self.nFN, 1)
    
    def __str__(self):
        return '%6d %6d %.2f %6d %.2f' % (self.nTP, self.nFP, self.FPRate(), self.nFN, self.FNRate())

def variantInTruth(variant, truth):
    if variant.getLoc() in truth:
        return truth[variant.getLoc()]
    else:
        return False

def sensitivitySpecificity(variants, truth):
    nTP, nFP = 0, 0    
    FPs = []
    for variant in variants:
        t = variantInTruth(variant, truth)
        if t:
            t.setField("FN", 0)
            variant.setField("TP", 1)
            nTP += 1
        else:
            nFP += 1
            #if variant.getPos() == 1520727:
            #    print "Variant is missing", variant
            FPs.append(variant)
    nFN = len(truth) - nTP
    return CallCmp(nTP, nFP, nFN), FPs


def compareCalls(calls, truthCalls):
    for variant in calls: variant.setField("TP", 0) # set the TP field to 0
    
    def compare1(name, cumulative):
        for left, right, selectedVariants in mapVariantBins(calls, 'QUAL', cumulative = cumulative):
            callComparison, theseFPs = sensitivitySpecificity(selectedVariants, truthCalls)
            #print selectedVariants[0]
            printFieldQual(name, left, right, selectedVariants, dephredScale(left), str(callComparison))
    
    print 'PER BIN nCalls=', len(calls)
    printFieldQualHeader("TP    FP  FPRate  FN  FNRate")
    compare1('TRUTH-PER-BIN', False)

    print 'CUMULATIVE nCalls=', len(calls)
    compare1('TRUTH-CUMULATIVE', True)
    
def randomSplit(l, pLeft):
    def keep(elt, p):
        if p < pLeft:
            return elt, None
        else:
            return None, elt
    data = [keep(elt, p) for elt, p in zip(l, map(lambda x: random.random(), l))]
    def get(i): return filter(lambda x: x <> None, [x[i] for x in data])
    return get(0), get(1)

def setup():
    global OPTIONS, header
    usage = "usage: %prog files.list [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--f", dest="fields",
                        type='string', default="QUAL",
                        help="Comma-separated list of fields (either in the VCF columns of as INFO keys) to use during optimization [default: %default]")
    parser.add_option("-t", "--truth", dest="truth",
                        type='string', default=None,
                        help="VCF formated truth file.  If provided, the script will compare the input calls with the truth calls.  It also emits calls tagged as TP and a separate file of FP calls")
    parser.add_option("", "--unFilteredTruth", dest="unFilteredTruth",
                        action='store_true', default=False,
                        help="If provided, the unfiltered truth calls will be used in comparisons [default: %default]")
    parser.add_option("-p", "--partitions", dest="partitions",
                        type='int', default=25,
                        help="Number of partitions to use for each feature.  Don't use so many that the number of variants per bin is very low. [default: %default]")
    parser.add_option("", "--maxRecordsForCovariates", dest="maxRecordsForCovariates",
                        type='int', default=200000,
                        help="Derive covariate information from up to this many VCF records.  For files with more than this number of records, the system downsamples the reads [default: %default]")
    parser.add_option("-m", "--minVariantsPerBin", dest="minVariantsPerBin",
                       type='int', default=10,
                        help="Don't include any covariates with fewer than this number of variants in the bin, if such a thing happens.  NEEDS TO BE FIXED")
    parser.add_option("-M", "--maxRecords", dest="maxRecords",
                       type='int', default=None,
                        help="Maximum number of input VCF records to process, if provided.  Default is all records")
    parser.add_option("-q", "--qMax", dest="maxQScore",
                        type='int', default=30,
                        help="The maximum Q score allowed for both a single covariate and the overall QUAL score [default: %default]")
    parser.add_option("-o", "--outputVCF", dest="outputVCF",
                        type='string', default=None,
                        help="If provided, a VCF file will be written out to this file [default: %default]")
    parser.add_option("", "--FNoutputVCF", dest="FNoutputVCF",
                        type='string', default=None,
                        help="If provided, VCF file will be written out to this file [default: %default]")
    parser.add_option("", "--titv", dest="titvTarget",
                        type='float', default=None,
                        help="If provided, we will optimize calls to the targeted ti/tv rather than that calculated from known calls [default: %default]")
    parser.add_option("-b", "--bootstrap", dest="bootStrap",
                       type='float', default=None,
                       help="If provided, the % of the calls used to generate the recalibration tables. [default: %default]")
    parser.add_option("-r", "--dontRecalibrate", dest="dontRecalibrate",
                        action='store_true', default=False,
                        help="If provided, we will not actually do anything to the calls, they will just be assessed [default: %default]")
   
    (OPTIONS, args) = parser.parse_args()
    if len(args) > 2:
        parser.error("incorrect number of arguments")
    return args

def assessCalls(file):
    print 'Counting records in VCF', file
    numberOfRecords = quickCountRecords(open(file))
    if OPTIONS.maxRecords <> None and OPTIONS.maxRecords < numberOfRecords:
        numberOfRecords = OPTIONS.maxRecords
    downsampleFraction = min(float(OPTIONS.maxRecordsForCovariates) / numberOfRecords, 1)
    header, allCalls = readVariants(file, OPTIONS.maxRecords, downsampleFraction=downsampleFraction)
    allCalls = list(allCalls)
    print 'Number of VCF records', numberOfRecords, ', max number of records for covariates is', OPTIONS.maxRecordsForCovariates, 'so keeping', downsampleFraction * 100, '% of records'
    print 'Number of selected VCF records', len(allCalls)
    
    titvtarget = OPTIONS.titvTarget
    if titvtarget == None:
        titvtarget = titv(selectVariants(allCalls, VCFRecord.isKnown))
    print 'Ti/Tv all  ', titv(allCalls)
    print 'Ti/Tv known', titv(selectVariants(allCalls, VCFRecord.isKnown))
    print 'Ti/Tv novel', titv(selectVariants(allCalls, VCFRecord.isNovel))
    
    return header, allCalls, titvtarget

def determineCovariates(allCalls, titvtarget, fields):
    if OPTIONS.bootStrap:
        callsToOptimize, recalEvalCalls = randomSplit(allCalls, OPTIONS.bootStrap)
    else:
        callsToOptimize = allCalls 

    recalOptCalls, covariates = optimizeCalls(callsToOptimize, fields, titvtarget)
    printCallQuals(list(recalOptCalls), titvtarget, 'OPTIMIZED CALLS')

    if OPTIONS.bootStrap:
        recalibatedEvalCalls = recalibrateCalls(recalEvalCalls, fields, covariates)
        printCallQuals(list(recalibatedEvalCalls), titvtarget, 'BOOTSTRAP EVAL CALLS')

    return covariates

def writeRecalibratedCalls(file, header, calls):
    if file:
        f = open(file, 'w')
        #print 'HEADER', header
        i = 0
        for line in formatVCF(header, calls):
            if i % 10000 == 0: print 'writing VCF record', i
            i += 1
            print >> f, line
        f.close()

def evaluateTruth(header, callVCF, truthVCF):
    print 'Reading truth file', truthVCF
    rawTruth = list(readVariants(truthVCF, maxRecords = None, decodeAll = False)[1])
    def keepVariant(t): 
        #print t.getPos(), t.getLoc()
        return OPTIONS.unFilteredTruth or t.passesFilters()
    truth = dict( [[v.getLoc(), v] for v in filter(keepVariant, rawTruth)])
    print len(rawTruth), len(truth)

    print 'Reading variants back in from', callVCF
    header, calls = readVariants(callVCF)
    calls = list(calls)
    
    print '--------------------------------------------------------------------------------'
    print 'Comparing calls to truth', truthVCF
    print ''

    compareCalls(calls, truth)

    writeRecalibratedCalls(callVCF, header, calls)

    if truth <> None and OPTIONS.FNoutputVCF:
        f = open(OPTIONS.FNoutputVCF, 'w')
        #print 'HEADER', header
        for line in formatVCF(header, filter( lambda x: not x.hasField("TP"), truth.itervalues())):
            print >> f, line
        f.close()

def main():
    args = setup()

    header, allCalls, titvTarget = assessCalls(args[0])
    if not OPTIONS.dontRecalibrate: 
        fields = OPTIONS.fields.split(',')
        covariates = determineCovariates(allCalls, titvTarget, fields)
        header, callsToRecalibate = readVariants(args[0], OPTIONS.maxRecords)
        RecalibratedCalls = recalibrateCalls(callsToRecalibate, fields, covariates)
        writeRecalibratedCalls(OPTIONS.outputVCF, header, RecalibratedCalls)
    else:
        printCallQuals(allCalls, titvTarget)
        OPTIONS.outputVCF = args[0]

    if len(args) > 1:
        evaluateTruth(header, OPTIONS.outputVCF, args[1])


PROFILE = False
if __name__ == "__main__":
    if PROFILE:
        import cProfile
        cProfile.run('main()', 'fooprof')
        import pstats
        p = pstats.Stats('fooprof')
        p.sort_stats('cumulative').print_stats(10)
        p.sort_stats('time').print_stats(10)
        p.sort_stats('time', 'cum').print_stats(.5, 'init')
    else:
        main()