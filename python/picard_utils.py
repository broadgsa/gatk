import farm_commands
import os.path
import sys
from optparse import OptionParser
import string
import re
import glob
import unittest
import itertools

#lanes = ["30JW3AAXX.6", "30KRNAAXX.1", "30KRNAAXX.6", "30PYMAAXX.5"]
#idsList = ['NA12843', 'NA19065', 'NA19064', 'NA18637']

lanes = ["30JW3AAXX.6", "30PYMAAXX.5", "30PNUAAXX.8", "30PPJAAXX.5"]
idsList = ['NA12843', 'NA18637', "NA19058", "NA12842"]
ids = dict(zip(lanes, idsList))
gatkPath = "~/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar"
ref = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"
analysis = "CombineDuplicates"

MERGE_BIN = '/seq/software/picard/current/bin/MergeSamFiles.jar'
SAMTOOLS_MERGE_BIN = '/seq/dirseq/samtools/current/samtools merge'
CALL_GENOTYPES_BIN = '/seq/software/picard/current/bin/CallGenotypes.jar'

def CollectDbSnpMatchesCmd(inputGeli, outputFile, lod): 
    return 'CollectDbSnpMatches.jar INPUT=%s OUTPUT=%s DBSNP=/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.dbsnp MINIMUM_LOD=%f' % (inputGeli, outputFile, lod)

def unique(l):
    return list(set(l))

def genotypes2heterozygosity(genotypes, nIndividuals = -1):
    def isHET(genotype):
        return genotype[0] <> genotype[1]
       
    if nIndividuals == -1:
        n = len(genotypes)
    else:
        n = nIndividuals
    
    hets = filter( isHET, genotypes )
    nhets = len(hets)
    # print genotypes, ' => hets', hets
    return [nhets / (1.0*n), nhets, n]

def genotypes2allelefrequencies(ref, genotypes, nIndividuals = -1):
    if nIndividuals == -1:
        n = len(genotypes)
    else:
        n = nIndividuals

    alleles = ''.join(genotypes)
    nChroms = 2 * n
    nCalledChroms = 2 * len(genotypes)
    nMissingChroms = nChroms - nCalledChroms
    nRefChroms = alleles.count(ref) + nMissingChroms
    nAltChroms = nChroms - nRefChroms
    p = float(nRefChroms) / nChroms
    q = float(nAltChroms) / nChroms

    #print 'genotypes', genotypes
    #print 'alleles', alleles
    #print 'nChroms', nChroms
    #print 'nCalledChroms', nCalledChroms
    #print 'nMissingChroms', nMissingChroms
    #print 'nRefChroms', nRefChroms
    #print 'nAltChroms', nAltChroms
    #print 'p, q', p, q
    
    assert p + q == 1
    
    return [p, q, n]

class PicardSNP:
    def __init__( self, loc, ref, polymorphism, heterozygosity, allelefrequencies, genotypes, sources):
        self.loc = loc
        self.ref = ref
        self.polymorphism = polymorphism
        self.heterozygosity = heterozygosity
        self.nIndividuals = allelefrequencies[2]
        self.allelefrequencies = allelefrequencies
        self.genotypes = genotypes
        self.sources = sources

    def refGenotype(self):
        return self.ref + self.ref

    def hetGenotype(self):
        return self.ref + self.alt()

    def homVarGenotype(self):
        return self.alt() + self.alt()

    def alt(self):
        return self.polymorphism[1]
    
    def het(self):
        return self.heterozygosity[0]
    
    def p(self):
        return self.allelefrequencies[0]
    
    def q(self):
        return self.allelefrequencies[1]
        
    def countGenotype(self, genotype):
        r = len(filter( lambda x: sorted(x) == sorted(genotype), self.genotypes ))
        #print 'countGenotype', genotype, self.genotypes, r
        return r

    def nRefGenotypes(self):
        return self.nIndividuals - self.nHetGenotypes() - self.nHomVarGenotypes()
        
    def nHetGenotypes(self):
        return self.countGenotype(self.hetGenotype()) 
    
    def nHomVarGenotypes(self):
        return self.countGenotype(self.homVarGenotype()) 

    def __str__(self):
        return '%s %s %s %s %s %s' % ( self.loc, self.ref, str(self.polymorphism), str(self.het()), str(self.allelefrequencies), str(self.genotypes))
        
def aggregatedGeliCalls2SNP( geliCallsAtSite, nIndividuals ):
    #print 'geliCallsAtSite', geliCallsAtSite
    loc = geliCallsAtSite[0]
    #print loc
    refBases = map( lambda call: call[2], geliCallsAtSite[1] )
    refBase = refBases[0]
    #print 'refBases', refBases

    genotypes = map( lambda call: ''.join(sorted(call[5])), geliCallsAtSite[1] )
    allBases = unique(''.join(genotypes) + refBase)
    #print 'All bases => ', allBases, genotypes 
    
    if len(allBases) > 2:
        print '*** WARNING, tri-state allele [ref=%s, all bases observed = %s] discovered at %s, ignoring the call' % ( refBase, ''.join(allBases), loc )
        return None

    #print 'genotypes', genotypes
    polymorphism = unique(list(refBase + genotypes[0]))
    if polymorphism[0] <> refBase: polymorphism.reverse()
    
    #print 'polymorphism', polymorphism
    genotype  = list(geliCallsAtSite[1][0][5])
    
    return PicardSNP(loc, refBase, polymorphism, genotypes2heterozygosity(genotypes, nIndividuals), genotypes2allelefrequencies(refBase, genotypes, nIndividuals), genotypes, [])

    #return '%s   %s %s 0.002747 -411.622578 -420.661738 0.000000 9.039160 364.000000 %d 1 0' % (loc, genotype[0], genotype[1], len(geliCallsAtSite))
    
def call2loc(call):
    return call[0] + ':' + call[1]

def aggregateGeliCalls( sortedGeliCalls ):
    #return [[loc, list(sharedCallsGroup)] for (loc, sharedCallsGroup) in itertools.groupby(sortedGeliCalls, call2loc)]
    return [[loc, list(sharedCallsGroup)] for (loc, sharedCallsGroup) in itertools.groupby(sortedGeliCalls, call2loc)]

def mergeBAMCmd( output_filename, inputFiles, mergeBin = MERGE_BIN, MSD = True, useSamtools = False, memLimit = '-Xmx4096m', compression_level = 1 ):
    if useSamtools:
        return SAMTOOLS_MERGE_BIN + ' ' + output_filename + ' ' + ' '.join(inputFiles)
    else:
        # use picard
        if type(inputFiles) <> list:
            inputFiles = list(inputFiles)
    
        MSDStr = ''
        if MSD: MSDStr = 'MSD=true'
    
        return 'java ' + memLimit + ' -jar ' + mergeBin + ' ' + MSDStr + ' AS=true COMPRESSION_LEVEL=' + str(compression_level) + ' SO=coordinate O=' + output_filename + ' VALIDATION_STRINGENCY=SILENT ' + ' I=' + (' I='.join(inputFiles))
        #return 'java -Xmx4096m -jar ' + mergeBin + ' AS=true SO=coordinate O=' + output_filename + ' VALIDATION_STRINGENCY=SILENT ' + ' I=' + (' I='.join(inputFiles))

def getPicardPath(lane, picardRoot = '/seq/picard/'):
    flowcell, laneNo = lane.split('.')
    filePat = os.path.join(picardRoot, flowcell, '*', laneNo, '*')
    dirs = glob.glob(filePat)
    print dirs
    if len(dirs) > 1:
        system.exit("Bad lane -- too many directories matching pattern " + filePat)
    return dirs[0]

def getReferenceGenotypeFileFromConcordanceFile(concordFile):
    # REFERENCE_GENOTYPES=/seq/references/reference_genotypes/hapmap/Homo_sapiens_assembly18/NA19058.geli
    p = re.compile('REFERENCE_GENOTYPES=([/.\w]+)')
    for line in open(concordFile):
        match = p.search(line)
        print 'Match is', line, match
        if match <> None:
            return match.group(1)
    return None

def hybridSelectionExtraArgsForCalling():
    return "TARGET_INTERVALS=/seq/references/HybSelOligos/thousand_genomes_alpha_redesign/thousand_genomes_alpha_redesign.targets.interval_list CALL_ZERO_COVERAGE_LOCI=true"

def callGenotypesCmd( inputBam, outputFilename, callGenotypesBin = CALL_GENOTYPES_BIN, options = ''):
    return "java -jar %s INPUT=%s OUTPUT=%s REFERENCE_SEQUENCE=%s CALLER_ALGORITHM=QUALITY_SCORE PRIOR_MODEL=SNP_FREQUENCY %s" % ( callGenotypesBin, inputBam, outputFilename, ref, options)

def concord(options, geli, output, genotypeFile):
    return ("java -jar /seq/software/picard/current/bin/CollectGenotypeConcordanceStatistics.jar OPTIONS_FILE=%s INPUT=%s OUTPUT=%s REFERENCE_GENOTYPES=%s MINIMUM_LOD=5.0" % ( options, geli, output, genotypeFile ) )

def readPicardConcordance(file):
    p = re.compile('HOMOZYGOUS_REFERENCE|HETEROZYGOUS|HOMOZYGOUS_NON_REFERENCE')
# CATEGORY    OBSERVATIONS        AGREE   DISAGREE        PCT_CONCORDANCE
# HOMOZYGOUS_REFERENCE    853     853     0       1
# HETEROZYGOUS    416     413     3  0.992788
# HOMOZYGOUS_NON_REFERENCE        235     231     4       0.982979
    types = [str, int, int, int, float]
    def parse1(line):
        return [f(x) for f, x in zip(types, line.split())]
    data = [parse1(line) for line in open(file) if p.match(line) <> None]
    return data

def splitPath(geli):
    root, filename = os.path.split(geli)
    s = filename.split('.')
    flowcellDotlane = '.'.join(s[0:2])
    ext = '.'.join(s[2:])
    return [root, flowcellDotlane, ext]

def read_dbsnp(dbsnp_matches):
    next = False
    for line in open(dbsnp_matches):
        s = line.split()
        if next:
            return s
        if len(s) > 0 and s[0] == "TOTAL_SNPS":
            next = True
    return []

# ------------------------------------------------------------------------------------------
#
# Unit testing!
#
# ------------------------------------------------------------------------------------------
class TestPicardUnils(unittest.TestCase):
    def setUp(self):
        import cStringIO
        dataString = """chr1    1105366 T       52      99      CT      10.559975       10.559975               -117.68 -93.107178      -116.616493     -45.536842      -88.591728      -92.043671      -20.964022      -116.014435     -44.473339      -31.523996
chr1    1105411 G       22      99      AG      12.484722       12.484722               -23.995817      -27.909206      -10.875731      -27.909206      -46.579994      -29.546518      -46.579994      -23.360453      -29.546518      -46.579994
chr1    1105411 G       29      99      AG      12.033216       12.033216               -30.641142      -34.376297      -14.483982      -35.457623      -53.197636      -32.525024      -53.498665      -26.517199      -33.606354      -54.579994
chr1    1105857 G       6       99      AG      7.442399        2.096584                -7.55462        -9.3608 -5.458036       -9.3608 -20.279993      -16.37723       -20.279993      -12.900434      -16.37723       -20.279993
chr1    1105857 G       7       99      AG      10.889011       1.795554                -7.406977       -9.514187       -5.611423       -9.514187       -23.879993      -19.97723       -23.879993      -16.500435      -19.97723       -23.879993
chr1    1106094 T       20      99      CT      6.747106        6.747106                -56.979992      -43.143734      -56.806652      -23.699236      -41.036522      -42.97039       -9.862975       -56.505623      -23.525892      -16.610081
chr1    1110294 G       42      99      AG      21.076984       21.076984               -44.285015      -49.702267      -17.869579      -49.831242      -80.276649      -48.442669      -80.404335      -38.946564      -48.571644      -80.405624
chr1    1111204 C       26      99      CT      11.364679       11.364679               -55.479992      -31.040928      -55.479992      -36.424099      -23.349712      -31.040928      -11.985033      -55.479992      -36.424099      -32.811741
chr1    1111204 C       29      99      TT      34.740601       3.890135                -52.282055      -48.646442      -52.704597      -17.954794      -44.565525      -48.464859      -13.715057      -52.100475      -17.773212      -9.824923
chr1    1111204 C       31      99      CT      18.784479       18.784479               -71.079994      -39.870823      -71.079994      -44.303268      -31.878578      -39.870823      -13.094099      -71.079994      -44.303268      -39.48679
"""
        dataFile = cStringIO.StringIO(dataString)
        self.nIndividuals = 10

        self.genotypesSets = aggregateGeliCalls(map( string.split, dataFile.readlines() ) )
        self.genotypes = map(lambda x: aggregatedGeliCalls2SNP(x, self.nIndividuals), self.genotypesSets )
        self.locs = ["chr1:1105366", "chr1:1105411", "chr1:1105857", "chr1:1106094", "chr1:1110294", "chr1:1111204"]
        self.nhets = [1, 2, 2, 1, 1, 2]
        self.altAlleles = [1, 2, 2, 1, 1, 4] 

        self.aaf = map( lambda x: (1.0*x) / (2 * self.nIndividuals), self.altAlleles )
        self.hets = map( lambda x: (1.0*x) / self.nIndividuals, self.nhets )

    def testGenotypesSize(self):
        self.assertEqual(len(self.genotypesSets), 6)

    def testGenotypes2Het(self):
        print 'testGenotypes2Het...'
        self.assertEqual(genotypes2heterozygosity(['AT']), [1, 1, 1])
        self.assertEqual(genotypes2heterozygosity(['AA']), [0, 0, 1])
        self.assertEqual(genotypes2heterozygosity(['TT']), [0, 0, 1])
        self.assertEqual(genotypes2heterozygosity(['AT', 'AT']), [1, 2, 2])
        self.assertEqual(genotypes2heterozygosity(['AA', 'AA']), [0, 0, 2])
        self.assertEqual(genotypes2heterozygosity(['AT', 'AA']), [0.5, 1, 2])
        self.assertEqual(genotypes2heterozygosity(['AT', 'TT']), [0.5, 1, 2])
        self.assertEqual(genotypes2heterozygosity(['AT', 'TT', 'AA']), [1.0/3, 1, 3])
        self.assertEqual(genotypes2heterozygosity(['AT', 'AT', 'AA']), [2.0/3, 2, 3])

        self.assertEqual(genotypes2heterozygosity(['AT', 'AT'], 10), [2.0/10, 2, 10])
        self.assertEqual(genotypes2heterozygosity(['AT', 'AA'], 10), [1.0/10, 1, 10])

    def testAlleleFreqs(self):
        print 'testAlleleFreqs...'
        self.assertEqual(genotypes2allelefrequencies('A', ['AT']), [0.5, 0.5, 1])
        self.assertEqual(genotypes2allelefrequencies('T', ['AT']), [0.5, 0.5, 1])
        self.assertEqual(genotypes2allelefrequencies('A', ['AA']), [1.0, 0.0, 1])
        self.assertEqual(genotypes2allelefrequencies('A', ['TT']), [0.0, 1.0, 1])
        
        self.assertEqual(genotypes2allelefrequencies('A', ['TT'], 2), [0.5, 0.5, 2])
        self.assertEqual(genotypes2allelefrequencies('A', ['AA'], 2), [1.0, 0.0, 2])
        self.assertEqual(genotypes2allelefrequencies('A', ['AT'], 2), [3.0/4, 1.0/4, 2])
        self.assertEqual(genotypes2allelefrequencies('A', ['AT'], 3), [5.0/6, 1.0/6, 3])
        self.assertEqual(genotypes2allelefrequencies('T', ['AT'], 3), [5.0/6, 1.0/6, 3])

        self.assertEqual(genotypes2allelefrequencies('A', ['AT', 'AT'], 3), [4.0/6, 2.0/6, 3])
        self.assertEqual(genotypes2allelefrequencies('A', ['AT', 'TT'], 3), [3.0/6, 3.0/6, 3])
        self.assertEqual(genotypes2allelefrequencies('A', ['AT', 'TT', 'AA']), [3.0/6, 3.0/6, 3])

    def testGenotypeSetLocs(self):
        for set, loc in zip(self.genotypesSets, self.locs):
            #print loc, set
            self.assertEqual(set[0], loc)

    def testGenotypeLocs(self):
        for genotype, loc in zip(self.genotypes, self.locs):
            self.assertEqual(genotype.loc, loc)

    def testGenotypeHets(self):
        print 'testGenotypeHets:'
        for genotype, het in zip(self.genotypes, self.hets):
            print ' => ', genotype, het
            self.assertEqual(genotype.het(), het)

    def testGenotypeAlleleFreqs(self):
        print 'testGenotypeAlleleFreqs:'
        for genotype, af in zip(self.genotypes, self.aaf):
            print ' => ', genotype, af
            self.assertEqual(genotype.allelefrequencies, [1 - af, af, self.nIndividuals])
            
    def testSplit(self):
        self.assertEqual(splitPath('/seq/picard/30GA9AAXX/C1-152_2008-10-23_2009-04-05/1/Solexa-8267/30GA9AAXX.1.observed_genotypes.geli'), ['/seq/picard/30GA9AAXX/C1-152_2008-10-23_2009-04-05/1/Solexa-8267', '30GA9AAXX.1', 'observed_genotypes.geli'])
        self.assertEqual(splitPath('/seq/picard/30GA9AAXX/C1-152_2008-10-23_2009-04-05/2/Solexa-8268/30GA9AAXX.2.observed_genotypes.geli'), ['/seq/picard/30GA9AAXX/C1-152_2008-10-23_2009-04-05/2/Solexa-8268', '30GA9AAXX.2', 'observed_genotypes.geli'])

if __name__ == '__main__':
    unittest.main()
