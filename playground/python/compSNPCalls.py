#!/usr/bin/env python

import sys, string
import os
import re
from itertools import *
from optparse import OptionParser
from memo import DiskMemoize, time_func

class ref_genome:
    """Reads reference genome in FASTA format into a dict"""

    def __init__(self, ref_genome_file):
        ref_genome.chr_offset = [[] for i in range(45)]
        chr_id = 0
        seq = ""
        for line in open(ref_genome_file):
            if line.startswith(">"):
                print line[1:],
                if line.startswith(">chrM"): # skip first > line
                    continue
                ref_genome.chr_offset[chr_id] = seq
                chr_id += 1
                seq = " " # make it 1 indexed instead of 0 indexed
                #if chr_id > 2:
                #    break
            else:
                seq += line.rstrip().upper()
        ref_genome.chr_offset[chr_id] = seq

    def __getitem__(self, key):
        return ref_genome.chr_offset[key]

AffyChr2Index = dict()
for i in range(1,23):
    AffyChr2Index[str(i)] = i
AffyChr2Index['MT'] = 0
AffyChr2Index['X'] = 23
AffyChr2Index['Y'] = 24

class GenotypeCall:
    #ref = time_func(ref_genome)("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
    def __init__( self, chr, pos, genotype, snpP, lod ):
        self.chr = chr
        self.pos = int(pos)
        self._site = chr + ':' + str(self.pos)
        self._isSNP = snpP
        self.genotype = string.join(map(string.upper, sorted(genotype)), '/')    # sorted list of bases at position
        self.lod = lod

    def refbase(self):
        return GenotypeCall.ref[AffyChr2Index[self.chr]][self.pos]

    def __hash__(self):
        return hash(self._site)

    def __eq__(self, other):
        return self._site == other._site

    def site(self): return self._site
    def isSNP(self): return self._isSNP
    
    def ref_het_hom(self):
        if self.genotype[0] <> self.genotype[2]:
            return 1 # het(erozygous non-ref)
        else:
            # homozygous something
            if self.genotype[0] == self.refbase:
                return 0 # ref
            else:
                return 2 # hom(ozygous non-ref)

    def isHET(self): return self.genotype[0] <> self.genotype[2]
    def isHOM(self): return self.genotype[0] == self.genotype[2]
            
    def __str__(self):
        return "%s:%s %s %s" % ( self.chr, self.pos, self.genotype, self.lod)
        
MAQGenotypeEncoding = {
    'A' : ['A', 'A'], 
    'C' : ['C', 'C'], 
    'T' : ['T', 'T'], 
    'G' : ['G', 'G'], 
    "M" : ['A', 'C'], 
    'K' : ['G', 'T'], 
    'Y' : ['C', 'T'], 
    'R' : ['A', 'G'], 
    'W' : ['A', 'T'], 
    'S' : ['C', 'G'], 
    'D' : ['A', 'G', 'T'], 
    'B' : ['C', 'G', 'T'], 
    'H' : ['A', 'C', 'T'], 
    'V' : ['A', 'C', 'G'], 
    'N' : ['A', 'C', 'G', 'T'] }

MAQ2STDChr = dict()
for i in range(1,23):
    MAQ2STDChr['chr'+str(i)] = str(i)
MAQ2STDChr['chrM'] = 'MT'
MAQ2STDChr['chrX'] = 'X'
MAQ2STDChr['chrY'] = 'Y'

def convertMAQChr(maqChr):
    #print 'convertMAQChr:', maqChr, MAQ2STDChr[maqChr]
    if maqChr in MAQ2STDChr:
        return MAQ2STDChr[maqChr]
    else:
        return '?'

def convertMAQGenotype( oneBaseCode ):
    return MAQGenotypeEncoding[oneBaseCode]

def internalReadSNPFile( parse1, filename ):
    result = []
    snps_extracted = 0
    for snp in imap( parse1, open(filename) ):
        if snp:
            result.append(snp)
            snps_extracted += 1
        if snps_extracted > OPTIONS.debug_lines:
            break

    print len(result),"genotypes extracted"
    return result

def snpMAP( snps ):
    #d = dict( map( lambda x: [x.site(), x], snps ) )
    d = dict()
    for snp in snps:
        d[snp.site()] = snp#d
    
    #print 'snps', snps, d
    return d

def overlappingSites( snps1, snps2 ):
    map1 = snpMAP(snps1)
    map2 = snpMAP(snps2)
    shared = set(map1.keys()) & set(map2.keys())
    print 'Number of snp1 records', len(map1)
    print 'Number of snp2 records', len(map2)
    print 'Number of shared sites', len(shared)
    print "\n".join(map(str,snps1))
    return shared

def readMAQSNPs(filename):
    # Each line consists of:
    #  chromosome
    #  position
    #  reference base
    #  consensus base
    #  Phred-like consensus quality
    #  read depth
    #  the average number of hits of reads covering this position
    #  the highest mapping quality of the reads covering the position
    #  the minimum consensus quality in the 3bp flanking regions at each side of the site (6bp in total)
    #  the second best call
    #  log likelihood ratio of the second best and the third best call
    #  and the third best call.
    # 
    # Also, note that:
    # 
    #      What do those "S", "M" and so on mean in the cns2snp output?
    # 		They are IUB codes for heterozygotes. Briefly:
    # 
    # 		M=A/C, K=G/T, Y=C/T, R=A/G, W=A/T, S=G/C, D=A/G/T, B=C/G/T, H=A/C/T, V=A/C/G, N=A/C/G/T
    def read1(line):
        formats = [str, int, str, str, int, int]
        vals = map( lambda f, x: f(x), formats, line.split()[0:6] )
        alignQual = vals[4]
        if alignQual >= (10*OPTIONS.lod):
            return GenotypeCall( convertMAQChr(vals[0]), vals[1], convertMAQGenotype(vals[3]), vals[2] <> vals[3], alignQual/10.0 )
        else:
            #print 'Filtering', alignQual, vals
            return False

    return internalReadSNPFile( read1, filename )

OPTIONS = None

def MerlinChr( index ):
    if index == 0: 
        return 'MT'
    elif index == 23: 
        return 'X'
    elif index == 24: 
        return 'Y'
    else:
        return str(index)
    
def readMerlinSNPs(filename):
    # 0:72 G GG 155.337967 0.000000 homozygous A:0 C:2 G:510 T:2 514 0 1 1 GG:-5.59 CG:-160.92 GT:-161.51 AG:-162.11 CT:-1293.61 CC:-1293.61 TT:-1294.19 AC:-1294.21 AT:-1294.80 AA:-1295.40 
    # 0:149 T CC 118.595886 1131.024696 homozygous-SNP A:2 C:442 G:1 T:7 452 0 1 1 CC:-24.21 CT:-142.81 AC:-156.33 CG:-156.96 TT:-1155.23 AT:-1159.41 GT:-1160.04 AA:-1173.26 AG:-1173.56 GG:-1174.20 
    # chr:pos ref genotype bestVsRef bestVsNextBest class ...
    def read1(line):
        formats = [lambda x: x.split(':'), str, sorted, float, float, str]
        vals = map( lambda f, x: f(x), formats, line.split()[0:6] )
        bestVsRef, bestVsNext = vals[3:5]
        isSNP = vals[5].find('-SNP') <> -1
        if bestVsRef >= OPTIONS.lod and isSNP:
            return GenotypeCall( MerlinChr(int(vals[0][0])), int(vals[0][1]) + 1, vals[2], isSNP, bestVsRef )
        else:
            return False

    return internalReadSNPFile( read1, filename )
    
def readSNPfile( filename, format ):
    formats = { 'merlin' : readMerlinSNPs, 'maq' : readMAQSNPs }
    if format.lower() in formats:
        return list(formats[format.lower()](filename))
    else:
        raise Exception('Unknown SNP file format ' + format)
    
def readAffyFile(filename):
    # chrom	position	genotype	probe_set_id	dbsnp_id
    # 1	84647761		TC	SNP_A-1780419	rs6576700
    # 5	156323558		GG	SNP_A-1780418	rs17054099
    def read1(line):
        formats = [str, int, sorted, str, str]
        vals = map( lambda f, x: f(x), formats, line.split() )
        
        try:
            chr = str(int(vals[0]))
        except:
            chr = convertMAQChr(vals[0])
        #print 'CHR', chr, vals[0]
        return GenotypeCall( chr, vals[1], vals[2], False, 100 )

    file = open(filename)
    file.readline()                                     # skip header
    #affyData = map( read1, file )
    affyData = []
    for index, line in enumerate(file):
        affyData.append(read1(line))
        if index > OPTIONS.debug_lines:
            break
        if index % 10000 == 0:
            print index
    # Give a chance to use list before creating dictionary
    return affyData

    #print "1111111"
    #return dict( zip( map( GenotypeCall.site, affyData ), affyData ) )

def equalSNPs( snp1, snp2 ):
    return snp1.genotype == snp2.genotype

# def concordance( truthSet, testSet, includeVector = None ):
#     # calculates a bunch of useful stats about the two 
#     # data genotype call sets above
#     
#     states = [[x,0] for x in ['tested', 'shared', 'shared-equal', 'test-snp', 'hom-ref', 'hom-snp', 'het-snp', 'het-ref']]
#     counts = dict(states)
#     
#     def incShared(state, equalP ):
#         counts[state] += 1
#         if equalP:
#             counts[state + '-equal'] += 1
#     
#     nTestSites = 0
#     for i, testSNP in izip( count(), testSet ):
#         if includeVector <> None and not includeVector[i]:
#             # we are skiping this site
#             continue
# 
#         nTestSites += 1
# 
#         if testSNP.isSNP():
#             counts['test-snp'] += 1
# 
#         #print testSNP.site()
#         if testSNP.site() in truthSet:
#             truth = truthSet[testSNP.site()]
#             eql = equalSNPs( testSNP, truth )
#             
#             incShared( 'shared', eql )
#             if testSNP.isSNP():
#                 if truth.isHOM(): incShared( 'hom-snp', eql )
#                 else: incShared( 'het-snp', eql )
#             else:
#                 if truth.isHOM(): incShared( 'hom-ref', eql )
#                 else: incShared( 'het-ref', eql )
#             
#         if OPTIONS.verbose and nTestSites % 100 == 0 and nSharedSites > 0:
#             #print nTestSites, nSharedSites, nEqualSites
#             print nTestSites, counts 
# 
#     counts['tested'] = nTestedSites
#             
#     return counts

class ConcordanceData:
    def __init__(self, name, file1count, file2count):
        self.name = name
        self.nFile1Sites = file1count    # num sites in file 1
        self.nFile2Sites = file2count    # num sites in file 1
        self.nSharedSites = 0   # num SNP pairs that map to same position on the genome
        self.nEqualSites = 0    # num SNPs pars with the same genotype

    def inc( self, truthSNP, testSNP ):
        self.nSharedSites += 1
        if equalSNPs( testSNP, truthSNP ): # if the genotypes are equal
            self.nEqualSites += 1
            
    def rate(self):
        return (100.0 * self.nEqualSites) / max(self.nSharedSites,1)        
        
    def __str__(self):
        return '%d %d %.2f' % ( self.nSharedSites, self.nEqualSites, self.rate() )

def concordance( truthSet, testSet, sharedSites = None ):
    # calculates a bunch of useful stats about the two 
    # data genotype call sets above
    #
    # The 2 calls in main work like this:
    # affy, snp1, snp1_snp2_shared
    # affy, snp2, snp1_snp2_shared
    
    nTestSites = 0

    # Now for each of the calls to concordance, we generate 3 sets:
    # - allData: all SNP1 sites that are also in Affy 
    allData = ConcordanceData('all', len(truthSet), len(testSet))
    # - sharedData: SNP1 sites that are also SNP2 sites that are alse in Affy
    sharedData = ConcordanceData('shared', len(truthSet), len(testSet))
    # - uniqueData: SNP1 sites that are not SNP2 sites but that are in Affy
    uniqueData = ConcordanceData('unique', len(truthSet), len(testSet))
    for i, testSNP in izip( count(), testSet ):
        nTestSites += 1
        if testSNP.site() in truthSet:
            truthSNP = truthSet[testSNP.site()]

            allData.inc( truthSNP, testSNP )
            if sharedSites <> None:
                if testSNP.site() in sharedSites:
                    sharedData.inc( truthSNP, testSNP )
                else:
                    uniqueData.inc( truthSNP, testSNP )

        if OPTIONS.verbose and nTestSites % 100000 == 0:
            #print nTestSites, nSharedSites, nEqualSites
            print nTestSites, allData, sharedData, uniqueData

    return nTestSites, allData, sharedData, uniqueData

# def concordance( truthSet, testSet, includeVector = None ):
#     # calculates a bunch of useful stats about the two 
#     # data genotype call sets above
#     
#     states = [[x,0] for x in ['tested', 'shared', 'test-snp', 'shared-hom-ref', 'shared-het-snp', 'shared-hom-snp']]
#     counts = dict(states)
#     
#     nTestSites = 0
#     nSharedSites = 0
#     nEqualSites = 0
#     for i, testSNP in izip( count(), testSet ):
#         nTestSites += 1
#         #print testSNP.site()
#         if testSNP.site() in truthSet:
#             nSharedSites += 1
#             if equalSNPs( testSNP, truthSet[testSNP.site()] ):
#                 nEqualSites += 1
#             #else:
#             #    print '~', testSNP, truthSet[testSNP.site()] 
#         if OPTIONS.verbose and nTestSites % 100000 == 0 and nSharedSites > 0:
#             #print nTestSites, nSharedSites, nEqualSites
#             print nTestSites, nSharedSites, nEqualSites, (100.0 * nEqualSites) / nSharedSites
#             
#     return [nTestSites, nSharedSites, nEqualSites, (100.0 * nEqualSites) / nSharedSites]


def printConcordanceResults( filename1, filename2, results, hasSharedSites = False ):
    nTestSites, allData, sharedData, uniqueData = results

    def print1(data):
        print '------------------------------------------------------------'
        print 'Concordance results', data.name, 'sites'
        print '  -> Genotype file1', filename1
        print '  -> Genotype file2', filename2
        print '  -> Number of tested sites (%s %s): %d' % (filename2, data.name, nTestSites)
        print '  -> Number of sites shared between files (%s %s): %d' % (filename2, data.name, data.nSharedSites)
        print '  ->   Percent sites in file1 shared (%s %s): %.2f' % (filename1, data.name, data.nSharedSites / (0.01 * data.nFile1Sites))
        print '  ->   Percent sites in file2 shared (%s %s): %.2f' % (filename1, data.name, data.nSharedSites / (0.01 * data.nFile2Sites))
        print '  -> Number of genotypically equivalent sites (%s %s): %d' % (filename2, data.name, data.nEqualSites)
        print '  -> Concordance rate (%s %s): %.2f' % (filename2, data.name, data.rate())
    
    print1(allData)
    if hasSharedSites:
        print1( sharedData ) # shared between SNP1 and SNP2
        print1( uniqueData ) # unique to SNP1 or to SNP2 only

def dump_shared_snps( affys, snp_list1, snp_list2 ):
    print len(affys), len(snp_list1), len(snp_list2)
    snps1 = snpMAP(snp_list1)
    snps2 = snpMAP(snp_list2)
    snp1_sites = set(snps1.keys()); print "SNP1s:",len(snp1_sites)
    snp2_sites = set(snps2.keys()); print "SNP2s:",len(snp2_sites)
    affy_sites = set(affys.keys()); print "Affys:",len(affy_sites)
    snp1or2_affy_sites = (snp1_sites | snp2_sites) & affy_sites
    snp1and2_affy_sites = (snp1_sites & snp2_sites) & affy_sites
    print "SNP 1 or 2 and Affy: ",len(snp1or2_affy_sites)
    print "SNP 1 and 2 and Affy:",len(snp1and2_affy_sites)

    fsnp = open ("snp.tab","w")
    print >>fsnp, "site lod1 lod2 lod1v2 gen1 gen2 genaff inc1 inc2 inc12 lodm genm incm ref_het_hom refbase"
    for site in snp1and2_affy_sites:
        snp1 = snps1[site]
        snp2 = snps2[site]
        affy = affys[site]

        print >>fsnp, "%-11s %5.2f %5.2f" % (site, snp1.lod, snp2.lod),
        try:
            snp1div2 = snp1.lod / snp2.lod
        except ZeroDivisionError:
            snp1div2 = 1000
        print >>fsnp, "%5.2f" % snp1div2,
        print >>fsnp, snp1.genotype, snp2.genotype, affy.genotype,
        print >>fsnp, "%1d %1d %1d" % (not equalSNPs(snp1, affy), not equalSNPs(snp2,affy), not equalSNPs(snp1,snp2)),

        # Calculte meta_lod from the two lods
        if snp1.genotype == snp2.genotype:
            meta_lod = snp1.lod + snp2.lod
            meta_genotype = snp1.genotype
        else:
            if snp1.lod > snp2.lod:
                meta_lod = snp1.lod
                meta_genotype = snp1.genotype
            else:
                meta_lod = snp2.lod
                meta_genotype = snp2.genotype
        meta_inc = meta_genotype != affy.genotype
        print >>fsnp, "%5.2f %3s %1d" % (meta_lod, meta_genotype, meta_inc),
        print >>fsnp, affy.ref_het_hom(),
        print >>fsnp, affy.refbase()

def intersection_union_snps( affy, snps1, snps2 ):
    map1 = snpMAP(snps1)
    map2 = snpMAP(snps2)
    shared_nonaffy_sites = (set(map1.keys()) & set(map2.keys())).difference(affy.keys())
    nonaffy_shared = [(map1[site].lod, map2[site].lod) for site in shared_nonaffy_sites]
    shared_affy_sites = set(map1.keys()) & set(map2.keys()) & set(affy.keys())

    #shared = []
    #for site in shared_sites:
    #    shared.append((map1[site], map2[site]))
    #print "Shared:",len( shared )
    #shared_x = [s[0].lod for s in shared]
    #shared_y = [s[1].lod for s in shared]

    both_corr = []
    snp1_corr = []
    snp2_corr = []
    neither_corr = []
    # given two bools telling whether snp1 and snp2 are correct, 
    # return the correct object
    #
    #             snp2 incorrect, snp2 correct
    truth_list = [[neither_corr, snp2_corr], # snp1 incorrect
                  [snp1_corr, both_corr]]    # snp1 correct

    for site in shared_affy_sites:
        snp1_true = equalSNPs(map1[site], affy[site])
        snp2_true = equalSNPs(map2[site], affy[site])
        truth_list[snp1_true][snp2_true].append( (map1[site].lod, map2[site].lod) )

    print "Beginning plot..."
    import rpy2.robjects as robj
    robj.r('X11(width=15, height=15)')
    XY_MAX = 25
    plots = ((nonaffy_shared, "gray45"),(both_corr, "black"), (snp1_corr, "red"), (snp2_corr, "green"), (neither_corr, "blue"))
    plots = ((both_corr, "black"), (snp1_corr, "red"), (snp2_corr, "green"), (neither_corr, "blue"))
    robj.r.plot([0, XY_MAX], [0, XY_MAX], \
                xlab = os.path.splitext(OPTIONS.snp1)[0].capitalize()+" LOD", \
                ylab = os.path.splitext(OPTIONS.snp2)[0].capitalize()+" LOD", \
                main = "Shared SNP site LODs", \
                xlim = robj.FloatVector([0,XY_MAX]), \
                ylim = robj.FloatVector([0,min(XY_MAX,25)]), \
                pch = 19, \
                cex = 0.0)
    for xy, color in plots:
        print "X"
        robj.r.points([pt[0] for pt in xy], [pt[1] for pt in xy], col=color, pch=19, cex=.3)
        print "Color:",color
        print "Len:",len(xy)
        #print "\n".join(["%25s %25s" % pt for pt in xy])

    #print "\n".join(["%25s %25s" % shared_snp for shared_snp in shared])
    #robj.r.plot(shared_x, shared_y, xlab=OPTIONS.format1.capitalize()+" LOD", ylab=OPTIONS.format2.capitalize()+" LOD", main="Shared SNP site LODs", col="black", xlim=robj.FloatVector([0,XY_MAX]), ylim=robj.FloatVector([0,min(XY_MAX,25)]), pch=19, cex=0.3)
    raw_input("Press enter to continue")
    
    return

    ss1 = set(snps1)
    ss2 = set(snps2)
    print "Shared:",len( ss1.intersection(ss2) )
    print "Snp1 only:",len( ss1.difference(ss2) )
    print "Snp2 only:",len( ss2.difference(ss1) )
    print "Snp1 total:",len( ss1 )
    print "Snp2 total:",len( ss2 )

def count_het_sites(snp_list):
    hets = 0
    for snp in snp_list:
        if snp.isHET():
            hets += 1
    print hets,"hets,",len(snp_list),"total,",
    print "%.1f" % (float(hets)/len(snp_list)*100)

def main(argv):
    global OPTIONS, ROOT

    usage = "usage: %prog --truth affy-truth-file --snp1 snpfile1 --snp2 snpfile2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--truth", dest="affy",
                        type="string", default=None,
                        help="Affy truth file")
    parser.add_option("-1", "--snp1", dest="snp1",
                        type="string", default=None,
                        help="Affy truth file")
    parser.add_option("-2", "--snp2", dest="snp2",
                        type="string", default=None,
                        help="Affy truth file")
    parser.add_option("", "--f1", dest="format1",
                        type="string", default=None,
                        help="File type of snpfile1")
    parser.add_option("", "--f2", dest="format2",
                        type="string", default=None,
                        help="File type of snpfile2")
    parser.add_option("-l", "--lod", dest="lod",
                        type="float", default=5.0,
                        help="Minimum LOD of confident SNP call in Merlin or Q (10x) in MAQ")
    parser.add_option("-v", "--verbose", dest="verbose",
                        action='store_true', default=False,
                        help="Verbose output")
    parser.add_option("-d", "--debug_lines", dest="debug_lines",
                        type='float', default=sys.maxint,
                        help="Number of input data lines to process for debugging")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 0:
        parser.error("incorrect number of arguments")

    if OPTIONS.affy == None:
        parser.error("No affy data specified")

    if OPTIONS.format1 == None:
        parser.error("First format cannot be none")

    if OPTIONS.snp2 <> None and OPTIONS.format2 == None:
        parser.error("snp2 file given but format was not specified")
        
    # Load reference genome
    #ref = ref_genome("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
    #sys.exit()

    if 1:
        print "Reading Affy truth data..."
        #readAffyFile2 = DiskMemoize( readAffyFile, "readAffyFile", global_deps = ["OPTIONS.lod"] )
        readAffyFile2 = time_func(readAffyFile)
        affy_list = readAffyFile2( filename=OPTIONS.affy )
        #count_het_sites(affy_list)
        #sys.exit()
        affy = dict( zip( map( GenotypeCall.site, affy_list ), affy_list ) )
        print 'Read affy truth data:'
        print '  -> number of genotyped loci', len(affy)

    #readSNPfile2 = DiskMemoize( readSNPfile, "readSNPfile", global_deps = ["OPTIONS.lod"] )
    readSNPfile2 = time_func(readSNPfile)
    print "Reading SNPs 1 file..."
    snps1 = readSNPfile2( filename=OPTIONS.snp1, format=OPTIONS.format1 )
    if OPTIONS.snp2 <> None:
        print "Reading SNPs 2 file..."
        snps2 = readSNPfile2( filename=OPTIONS.snp2, format=OPTIONS.format2 )
       
        dump_shared_snps( affy, snps1, snps2 )
    #intersection_union_snps( affy, snps1, snps2 )
    #sys.exit()

    sharedSites = None
    if OPTIONS.snp2 <> None:
        sharedSites = overlappingSites( snps1, snps2 )
        results1 = concordance( affy, snps1, sharedSites )
        printConcordanceResults( OPTIONS.affy, OPTIONS.snp1, results1, True )
        results2 = concordance( affy, snps2, sharedSites )
        printConcordanceResults( OPTIONS.affy, OPTIONS.snp2, results2, True )

    else:
        results = concordance( affy, snps1 )
        printConcordanceResults( OPTIONS.affy, OPTIONS.snp1, results, sharedSites )    

if __name__ == "__main__":       
    main(sys.argv)
