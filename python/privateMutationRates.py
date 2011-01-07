import sys
from optparse import OptionParser
from itertools import *
import random

# a simple script that does:
# 1 -- generates a master set of variants following the neutral expectation from a single big population
# 2 -- randomly generates M individuals with variants and genotypes sampled as expected from the big population of variants
# 3 -- writes out the genotypes of these individuals, and their allele frequency
def main():
    global OPTIONS
    usage = "usage: %prog [options] outputFile"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-N", "", dest="bigPopSize",
                        type='int', default=1000,
                        help="")         

    parser.add_option("-M", "", dest="smallPopSize",
                        type='int', default=100,
                        help="")         

    parser.add_option("-K", "", dest="nHetsPerSample",
                        type='int', default=1000,
                        help="")         

    parser.add_option("", "--maxMAF", dest="maxMAF",
                        type='float', default=None,
                        help="")         
         
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("Takes no arguments")

    random.seed(10000)
    genotypes = simulateSeqExpt(OPTIONS.bigPopSize, OPTIONS.smallPopSize, OPTIONS.nHetsPerSample)
    printGenotypes(genotypes, open(args[0] + ".genotypes", 'w'))
    printAFS(genotypes, open(args[0] + ".afs", 'w'))

class Variant:
    def __init__(self, id, trueAC, trueAN):
        self.id = "%d.%d" % ( trueAC, id )
        self.trueAC = trueAC
        self.trueAN = trueAN

        q = self.af()
        p = 1 - q
        self.hw = [p * p, 2 * p * q, q * q]
        
    def __str__(self):
        return "[V %s ac=%d an=%d af=%.2f]" % (self.id, self.trueAC, self.trueAN, self.af())
    __repr__ = __str__
        
    def af(self):
        return self.trueAC / (1.0*self.trueAN)
        
    def hwe(self): # returns phomref, phet, phomvar
        return self.hw

def simulateSeqExpt(bigPopSize, smallPopSize, nHetsPerSample):
    """Master runner function"""
    trueAFS = makeAFS(bigPopSize, nHetsPerSample)
    
    variants = AFStoVariants(trueAFS, bigPopSize)
    
    # returns a list of variants per sample
    genotypes = genotypeSamples(variants, smallPopSize)
    
    return genotypes
    
def makeAFS(nSamples, nHetsPerSample):
    """Generates allele frequency spectrum counts for nsamples and nHetsPerSample from neutral expectation"""
    nTotalVariants = nHetsPerSample * sum([1 / (1.0*i) for i in range(1, nSamples * 2 + 1)])
    AFSCounts = [int(round(nHetsPerSample / (1.0*i))) for i in range(1, nSamples * 2 + 1)]
    print AFSCounts
    print nTotalVariants
    print sum(AFSCounts)
    return AFSCounts

def AFStoVariants(trueAFS, bigPopSize):
    """Converts an allele frequency spectrum to specific named Variant objects"""
    variants = []

    nChromosomes = 2 * bigPopSize 
    for ac in range(len(trueAFS)):
        af = (1.0*ac) / nChromosomes
        if OPTIONS.maxMAF == None or af <= OPTIONS.maxMAF: 
            for j in range(trueAFS[ac]):
                v = Variant(j, ac+1, nChromosomes)
                #print ac, j, v
                variants.append(v)
        else:
            print 'Skipping AC', ac, ' / ', nChromosomes, 'beyond max MAF', OPTIONS.maxMAF

    return variants

# returns a list of variants per sample
def genotypeSamples(variants, nSamples):
    """Given a list of variants, generates nSamples genotypes"""
    return [genotypeSample(samplei, variants) for samplei in range(nSamples)]

def genotypeSample(id, variants):
    """Generate a single set of genotypes for a single using the list of variants"""
    print 'Genotyping sample', id
    genotypes = []
    for v in variants:
        pHomRef, pHet, pHomVar = v.hwe()
        r = random.random()
        if r > pHomRef:             # are we not reference?
            if r > pHomRef + pHet:  # are we hom var?
                count = 2
            else:
                count = 1
            #print (r, v.af(), pHomRef, pHet, pHomVar, count)
            genotypes.append([v, count])

    return genotypes    

def printGenotypes(sampleGenotypes, out):
    print >> out, "\t".join(["sample", "id", "ac", "an", "g"])
    for sample, i in izip(sampleGenotypes, count(len(sampleGenotypes))):
        for v, g in sample:
            print >> out, "\t".join(map(str, [i-1, v.id, v.trueAC, v.trueAN, g]))

def printAFS(sampleGenotypes, out):
    print >> out, "\t".join(["id", "true.ac", "true.an", "true.af", "small.ac", "small.an", "small.af"])
    counts = dict()

    smallAN = len(sampleGenotypes) * 2
    for sample in sampleGenotypes:
        for v, g in sample:
            if v not in counts: counts[v] = 0
            counts[v] = counts[v] + g
            

    for v, smallAC in counts.iteritems():
        print >> out, "\t".join(map(str, [v.id, v.trueAC, v.trueAN, v.af(), smallAC, smallAN, smallAC / (1.0*smallAN)]))


if __name__ == "__main__":
    main()
