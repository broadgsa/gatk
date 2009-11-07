import os.path
import sys
from optparse import OptionParser
from vcfReader import *
#import pylab
from itertools import *
import math


class RecalibratedCall:
    def __init__(self, call, features):
        self.call = call
        self.features = dict([[feature, None] for feature in features])
        
    def recalFeature( self, feature, FPRate ):
        assert self.features[feature] == None # not reassigning values
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
        jointTPRate = math.pow(10, logJointTPRate)
        #print logTPRates
        #print logJointTPRate, jointTPRate
        return 1 - jointTPRate
        
    def featureStringList(self):
        return ','.join(map(lambda feature: '%s=Q%d' % (feature, self.getFeature(feature, '*', True)), self.features.iterkeys()))        
        
    def __str__(self):
        return '[%s: %s => Q%d]' % (str(self.call), self.featureStringList(), phredScale(self.jointFPErrorRate()))

def readVariants( file ):
    counter = OPTIONS.skip

    def parseVariant(args):
        VCF, counter = args
        if counter % OPTIONS.skip == 0:
            return VCF
        else:
            return None

    return filter(None, map(parseVariant, lines2VCF(open(file))))

def selectVariants( variants, selector = None ):
    if selector <> None:
        return filter(selector, variants)
    else:
        return variants

def titv(variants):
    ti = len(filter(VCFRecord.isTransition, variants))
    tv = len(variants) - ti
    titv = ti / (1.0*max(tv,1))

    return titv, ti, tv


def gaussian(x, mu, sigma):    
    constant = 1 / math.sqrt(2 * math.pi * sigma**2)
    exponent = -1 * ( x - mu )**2 / (2 * sigma**2)
    return constant * math.exp(exponent)

DEBUG = False

# if target = T, and FP calls have ti/tv = 0.5, we want to know how many FP calls
# there are in N calls with ti/tv of X.  
# 
def titvFPRateEstimate(variants, target):
    titvRatio, ti, tv = titv(variants)
    
    # f <- function(To,T) { (To - T) / (1/2 - T) + 0.001 }
    def theoreticalCalc():
        if titvRatio >= target:
            FPRate = 0
        else:
            FPRate = (titvRatio - target) / (0.5 - target)
        FPRate = min(max(FPRate, 0), 1)
        TPRate = max(min(1 - FPRate, 1 - dephredScale(OPTIONS.maxQScore)), dephredScale(OPTIONS.maxQScore))
        print 'FPRate', FPRate, titvRatio, target
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
    if nVariants > 0:
        impliedNoErrors = nVariants * FPRate
        calcTiTv = (impliedNoErrors * 0.5 + target * (nVariants-impliedNoErrors)) / nVariants
    else:
        impliedNoErrors, calcTiTv = 0, 0
    
    if DEBUG: print ':::', nVariants, titvRatio, target, ti, tv, FPRate, impliedNoErrors, calcTiTv
    
    return titvRatio, FPRate, impliedNoErrors
    
def phredScale(errorRate):
    return -10 * math.log10(max(errorRate, 1e-10))

def dephredScale(qscore):
    return math.pow(10, qscore / -10)

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

def calculateBins(variants, field, minValue, maxValue, rangeValue, partitions):
    sortedVariants = sorted(variants, key = lambda x: x.getField(field))
    sortedValues = map(lambda x: x.getField(field), sortedVariants)
    
    targetBinSize = len(variants) / (1.0*partitions)
    print 'targetBinSize', targetBinSize
    uniqBins = groupby(sortedValues)
    binsAndSizes = map(lambda x: [x[0], len(list(x[1]))], uniqBins)
    #print binsAndSizes

    def bin2Break(bin): return [bin[0], bin[0], bin[1]]
    bins = [bin2Break(binsAndSizes[0])]
    for bin in binsAndSizes[1:]:
        #print 'Breaks', bins
        curSize = bin[1]
        prevSize = bins[-1][2]
        if curSize + prevSize > targetBinSize:
            bins.append(bin2Break(bin))
        else:
            bins[-1][1] = bin[0]
            bins[-1][2] += curSize

    return bins


def calculateBinsLinear(variants, minValue, maxValue, rangeValue, partitions):
    breaks = list(frange6(minValue, maxValue, rangeValue / partitions))
    if breaks[len(breaks)-1] <> maxValue:
        breaks = breaks + ['*']
    return zip(breaks, map( lambda x: x - 0.001, breaks[1:]))

def fieldRange(variants, field):
    values = map(lambda v: v.getField(field), variants)
    minValue = min(values)
    maxValue = max(values)
    rangeValue = maxValue - minValue
    bins = calculateBins(variants, field, minValue, maxValue, rangeValue, OPTIONS.partitions)
    return minValue, maxValue, rangeValue, bins

def printFieldQual( left, right, variants, titv, FPRate, nErrors ):
    print '  %s nVariants=%8d titv=%.2f FPRate=%.2e Q%d' % (binString(left, right), len(variants), titv, FPRate, phredScale(FPRate))

def binString(left, right):
    leftStr = str(left)
    if type(left) == float: leftStr = "%.2f" % left
    rightStr = "%5s" % str(right)
    if type(right) == float: rightStr = "%.2f" % right
    return '%8s - %8s' % (leftStr, rightStr)
    
#
#
#
def optimizeCalls(variants, fields, titvTarget):
    recalCalls = calibrateFeatures(variants, fields, titvTarget)

    newCalls = list()
    for recalCall in recalCalls.itervalues():
        originalQual = recalCall.call.getField('QUAL') 
        recalCall.call.setField('QUAL', int(round(phredScale(recalCall.jointFPErrorRate())))) 
        recalCall.call.setField('OQ', originalQual)
        newCalls.append(recalCall.call)

    for recalCall in islice(recalCalls.itervalues(), 10):
        print recalCall
    
    print '--------------------------------------------------------------------------------'
    print 'RECALIBRATED CALLS'
    #newCalls = [x.call for x in recalCalls.itervalues()]
    calibrateFeatures(newCalls, ['QUAL'], titvTarget, updateCalls = False, printCall = True )

    return newCalls

def all( p, l ):
    for elt in l:
        if not p(elt): return False
    return True

def variantBinsForField(variants, field):
    if not all( lambda x: x.hasField(field), variants):
        raise Exception('Unknown field ' + field)
    
    minValue, maxValue, range, bins = fieldRange(variants, field)
    print 'Field range', minValue, maxValue, range
    print 'Partitions', bins
    return bins

def mapVariantBins(variants, field):
    bins = variantBinsForField(variants, field)
    
    def variantsInBin(bin):
        left, right = bin[0:2] 
        def select( variant ): return variant.getField(field) >= left and (right == '*' or variant.getField(field) <= right)
        return left, right, selectVariants(variants, select)
        
    return imap( variantsInBin, bins )

def calibrateFeatures(variants, fields, titvTarget, updateCalls = True, printCall = False):
    if updateCalls: recalCalls = dict([[variant, RecalibratedCall(variant, fields)] for variant in variants]) 
    
    for field in fields:
        print 'Optimizing field', field
        
        titv, FPRate, nErrors = titvFPRateEstimate(variants, titvTarget)
        print 'Overall FRRate:', FPRate, nErrors, phredScale(FPRate)

        for left, right, selectedVariants in mapVariantBins(variants, field):
            if len(selectedVariants) > max(OPTIONS.minVariantsPerBin,1):
                titv, FPRate, nErrors = titvFPRateEstimate(selectedVariants, titvTarget)
                if updateCalls: 
                    for variant in selectedVariants: recalCalls[variant].recalFeature(field, FPRate)
                if printCall:
                    for call in selectedVariants:
                        if titv < 0.5:
                            print call

                printFieldQual( left, right, selectedVariants, titv, FPRate, nErrors )

    if updateCalls: 
        return recalCalls
    else:
        return None

class CallCmp:
    def __init__(self, nTP, nFP, nFN):
        self.nTP = nTP
        self.nFP = nFP
        self.nFN = nFN
    
    def FPRate(self):
        return (1.0*self.nFP) / max(self.nTP + self.nFP, 1)
    
    def __str__(self):
        return 'TP=%6d FP=%6d FPRate=%.2f FN=%6d' % (self.nTP, self.nFP, self.FPRate(), self.nFN)

def variantInTruth(variant, truth):
    return variant.getLoc() in truth

def sensitivitySpecificity(variants, truth):
    nTP, nFP = 0, 0    
    for variant in variants:
        if variantInTruth(variant, truth):
            nTP += 1
        else:
            if OPTIONS.printFP: print 'FP:', variant
            nFP += 1
    nFN = len(truth) - nTP
    return CallCmp(nTP, nFP, nFN)


def compareCalls(optimizedCalls, truthCalls):
    for left, right, selectedVariants in mapVariantBins(optimizedCalls, 'QUAL'):
        callComparison = sensitivitySpecificity(selectedVariants, truthCalls)
        print binString(left, right), 'titv=%.2f' % titv(selectedVariants)[0], callComparison
    


def main():
    global OPTIONS
    usage = "usage: %prog files.list [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--f", dest="fields",
                        type='string', default="QUAL",
                        help="Comma-separated list of fields to exact")
    parser.add_option("-t", "--truth", dest="truth",
                        type='string', default=None,
                        help="VCF formated truth file")
    parser.add_option("-p", "--partitions", dest="partitions",
                        type='int', default=25,
                        help="Number of partitions to examine")
    parser.add_option("-s", "--s", dest="skip",
                        type='int', default=1,
                        help="Only work with every 1 / skip records")
    parser.add_option("-m", "--minVariantsPerBin", dest="minVariantsPerBin",
                       type='int', default=10,
                        help="")
    parser.add_option("-q", "--qMax", dest="maxQScore",
                        type='int', default=30,
                        help="")
    parser.add_option("", "--titv", dest="titvTarget",
                        type='float', default=None,
                        help="If provided, we will optimize calls to the targeted ti/tv rather than that calculated from known calls")
    parser.add_option("", "--fp", dest="printFP",
                       action='store_true', default=False,
                        help="")
         
    (OPTIONS, args) = parser.parse_args()
    if len(args) > 2:
        parser.error("incorrect number of arguments")

    fields = OPTIONS.fields.split(',')
    calls = readVariants(args[0])
    print 'Read', len(calls), 'calls'
    
    if OPTIONS.titvTarget == None:
        OPTIONS.titvTarget = titv(calls, VCFRecord.isKnown)
    print 'Ti/Tv all  ', titv(calls)
    print 'Ti/Tv known', titv(selectVariants(calls, VCFRecord.isKnown))
    print 'Ti/Tv novel', titv(selectVariants(calls, VCFRecord.isNovel))

    optimizedCalls = optimizeCalls(calls, OPTIONS.fields.split(","), OPTIONS.titvTarget)

    if len(args) > 1:
        truthFile = args[1]
        print 'Reading truth file', truthFile
        truth = dict( [[v.getLoc(), v] for v in readVariants(truthFile)])
        compareCalls(optimizedCalls, truth)

if __name__ == "__main__":
    main()