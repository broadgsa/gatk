#!/util/bin/python

from optparse import OptionParser
import os
import sys
#import set
import random
import string

from SAM import *

def all(iterable):
    for element in iterable:
        if not element:
            return False
    return True

def any(iterable):
    for element in iterable:
        if element:
            return True
    return False

OPTIONS = None

bases = set(['a', 't', 'g', 'c'])

class Mutation:
    """Representation of a change from one sequence to another""" 

    SNP = "SNP"
    INSERTION = "Insertion"
    DELETION = "Deletion"
    TYPES = set([SNP, INSERTION, DELETION])
    
    CIGARChars = { SNP : 'M', INSERTION : 'I', DELETION : 'D' }

    def __init__( self, contig, pos, type, original, mutated ):
        # example: chr20 123 SNP "a" "t"        <- a to t at base 123 on chr20
        # example: chr20 123 INSERTION "" "atg" <- insertion of 3 bases starting at pos 123
        # example: chr20 123 DELETION  "atg" "" <- deletion of 3 bases starting at pos 123
        self.contig = contig    # contig of the mutation
        self._pos = pos          # position of the change in contig coordinates
        self._type = type        # one of the TYPES
        assert(self._type in Mutation.TYPES)
        self.original = original    # String of the original genotype
        self.mutated  = mutated     # String of the mutated genotype

    def pos(self, offset=0): return self._pos - offset
    def type(self): return self._type
    def __len__(self):
        return max( len(self.mutated), len(self.original) )

    def mutString( self ):
        if self._type == Mutation.SNP:
            return 'P:%s->%s' % (self.original, self.mutated)
        elif self._type == Mutation.INSERTION:
            return 'I:%s' % (self.mutated)
        elif self._type == Mutation.DELETION:
            return 'D:%s' % (self.original)
        else:
            raise Exception('Unexpected mutation type ' + self._type)

    def CIGARChar( self ):
        return Mutation.CIGARChars[self.type()]
            
    def __str__( self ):
        return '%s:%d %s' % (self.contig, self._pos, self.mutString())
        
    def __repr__( self ):
        return 'Mutation(%s, %d, %s)' % (self.contig, self._pos, self.mutString())

    def execute( self, offset, refseq, mutseq ):
        p = self.pos(offset)
        if self.type() == Mutation.SNP:
            mutseq[p] = self.mutated
        elif self.type() == Mutation.INSERTION:
            refseq[p:p] = string.join(['-'] * len(self),'')
            mutseq[p:p] = self.mutated
        elif self.type() == Mutation.DELETION:
            mutseq[p:p+len(self)] = string.join(['-'] * len(self), '')
            refseq.extend(self.mutated)
            mutseq.extend(self.mutated)
            
        #print refseq, mutseq
        return refseq, mutseq
        
# ---------------------------------------------------------------------------
#
# Code for dealing with CIGAR strings
#
# ---------------------------------------------------------------------------
def flatten(l):
    out = []
    for item in l:
        if isinstance(item, (list, tuple)):
            out.extend(flatten(item))
        else:
            out.append(item)
    return out

def mutationsCIGAR( readStart, alignStart, readLen, mutations ):
    if len(mutations) == 1 and mutations[0].type() <> Mutation.SNP:
        mut = mutations[0]
        internalPos = mut.pos() - alignStart

        leftLen = max(internalPos - readStart,0)
        mutLen = min(internalPos + len(mut) - readStart, readLen - leftLen, len(mut))
        rightLen = readLen - leftLen - mutLen
        #print mut, readStart, alignStart, internalPos, leftLen, mutLen, rightLen
        #print
        
        l = [ (leftLen, 'M'), (mutLen, mut.CIGARChar()), (rightLen, 'M') ]
        l = filter( lambda x: x[0] > 0, l )
        return string.join(map(str, flatten(l)), '')
        
    if not all( map( lambda mut: mut.type() == Mutation.SNP, mutations ) ):
        # bail
        raise Exception('Currently only supports multiple SNPs CIGARS')

    return str(readLen) + 'M' 

# def mutationsCIGARBAD( refStart, readStart, readStop, mutations ):
#     sortedMutations = sorted(mutations, key=Mutation.pos)
#     CIGAR = []
#     pos = readStart
#     
#     for mutation in sortedMutations:
#         internalMutationStart = mutation.pos() - refStart
#         if internalMutationStart >= readStart and internalMutationStart < readStop
#             delta = mutation.pos() - pos
#             if not OPTIONS.quiet:
#                 print 'mutationsCIGAR', pos, CIGAR, delta, mutation
#             if mutation.type() <> Mutation.SNP and delta > 0:
#                 CIGAR.append([delta, 'M'])
#                 pos = mutation.pos()
#             
#             if mutation.type == Mutation.INSERTION:
#                 CIGAR.append(['I', len(mutation)])
#             if mutation.type == Mutation.DELETION:
#                 CIGAR.append(['D', len(mutation)])
# 
#     delta = refSeqPos + refLen - pos
#     if not OPTIONS.quiet:
#         print 'mutationsCIGAR', pos, CIGAR, delta, mutation
#     if delta > 0:
#         CIGAR.append([delta, 'M'])
# 
#     s = string.join(map( str, flatten(CIGAR) ), '')
#     print 'CIGAR:', CIGAR, s
#     return s

# ---------------------------------------------------------------------------
#
# mutating the reference
#
# ---------------------------------------------------------------------------
def executeMutations( mutations, seqIndex, refseq, mutseq ):
    for mutation in mutations:
        internalSite = mutation.pos() - seqIndex
        refseq, mutseq = mutation.execute( seqIndex, refseq, mutseq )
    return refseq, mutseq
    
def isBadSequenceForSim(seq):
    if any( map( lambda b: b not in bases, seq.tostring().lower() ) ):
        return True
    else:
        return False
    

def mutateReference(ref, contig, mutSite, mutType, mutParams, nBasesToPad):
    #print 'mutateReference'

    # start and stop    
    refStartIndex = mutSite - nBasesToPad + 1
    internalStartIndex = nBasesToPad - 1

    if OPTIONS.mutationType == 'SNP' or OPTIONS.mutationType == 'NONE':
        refStopIndex = mutSite + nBasesToPad
    else:
        refStopIndex = mutSite + nBasesToPad - 1
    
    refsub = ref[refStartIndex : refStopIndex]
    mutseq = refsub.tomutable()
    refsub = refsub.tomutable()
    mutations = []

    def success():
        return True, mutations, refsub, mutseq
    def fail():
        return False, [], None, None 

    if refStartIndex < 0: 
        return fail()
    if isBadSequenceForSim(refsub):
        #print 'rejecting seq', refsub.tostring()
        #print map( lambda b: b not in bases, refsub.tostring().lower() )
        return fail()

    if not OPTIONS.quiet:
        print ref
        print refsub

    otherSites = set(range(refStartIndex, refStopIndex)) - set([mutSite])
    # print 'otherSites', otherSites
    for site in [mutSite] + random.sample(otherSites, max(OPTIONS.mutationDensity-1,0)):
        internalSite = site - refStartIndex
        if mutType == 'NONE' or OPTIONS.mutationDensity < 1:
            pass
        elif mutType == 'SNP':
            mutBase = ref[site].lower()
            otherBases = bases - set(mutBase)
            otherBase = random.choice(list(otherBases))
            assert mutBase <> otherBase
            if not OPTIONS.quiet:
                print mutBase, ' => ', otherBase
            mutations.append( Mutation( contig, site, Mutation.SNP, mutBase, otherBase ) )
        elif mutType == 'INSERTION':
            inSeq = string.join([ random.choice(list(bases)) for i in range( OPTIONS.indelSize ) ], '')
            mutation = Mutation( contig, site, Mutation.INSERTION, '', inSeq )
            mutations.append( mutation )
        elif mutType == 'DELETION':
            inSeq = string.join([ random.choice(list(bases)) for i in range( OPTIONS.indelSize ) ], '')
            mutation = Mutation( contig, site, Mutation.DELETION, '', ref[refStopIndex:refStopIndex+OPTIONS.indelSize])
            mutations.append( mutation )
        else:
            raise Exception('Unexpected mutation type', mutType)

    # process the mutations
    refsub, mutseq = executeMutations( mutations, refStartIndex, refsub, mutseq )

    return success()

# version 1.0 -- prototype
# def mutateReference(ref, mutSite, mutType, mutParams, nBasesToPad):
#     #print 'mutateReference'
#     refStartIndex = mutSite - nBasesToPad + 1
#     refStopIndex = mutSite + nBasesToPad
#     internalStartIndex = nBasesToPad - 1
# 
#     if refStartIndex < 0:
#         return False, None, None
# 
#     refsub = ref[refStartIndex : refStopIndex]
#     print ref
#     print refsub
# 
#     if mutType == 'SNP' or mutType == 'None':
#         #print '  In IF'
#         mutBase = ref[mutSite]
#         bases = set(['A', 'T', 'G', 'C'])
#         if mutBase in bases:
#             otherBases = bases - set(mutBase)
#             otherBase = random.choice(list(otherBases))
#             assert mutBase <> otherBase
#             mutseq = refsub.tomutable()
# 
#             if mutType == 'SNP':
#                 print mutBase, ' => ', otherBase
#                 mutseq[internalStartIndex] = otherBase
# 
#             print refsub
#             print mutseq
#             align = Bio.Align.Generic.Alignment(Gapped(IUPAC.unambiguous_dna, '-'))
#             align.add_sequence("ref", refsub.tostring())
#             align.add_sequence("mut", mutseq.tostring())
#             print str(align)
# 
#             return True, refsub, mutseq
# 
#     return False, None, None
            
# ---------------------------------------------------------------------------
#
# sample reads from alignments
#
# ---------------------------------------------------------------------------
def accumulateBases( mutSeq, start, readLen ):
    count = 0
    for stop in range(start, len(mutSeq)):
        if notDash(mutSeq[stop]):
            count += 1
        #print stop, count
        if count == readLen:
            break
            
    stop += 1
    read = string.join(filter( notDash, mutSeq[start:stop] ), '')

    return read, stop

# doesn't support paired reads
#
# def sampleReadsFromAlignment(refSeq, mutSeq, alignStart, readLen, nReads, mutations):
#     #print 'REF:', refSeq.tostring()
#     #print 'MUT:', mutSeq.tostring()
#     #print ''
#     
#     # we are assuming that mutSeq is exactly 1 + 2 * readLen in size, which the mutation
#     # right in the middle
#     lastGoodStartSite = readLen - 1
#     if not OPTIONS.quiet:
#         print refSeq.tostring(), 'ref'
#         print mutSeq.tostring(), 'mut'
#     
#     # we can potentially start at any site in mutSeq
#     nPotentialStarts = len(mutSeq) - readLen + 1
#     potentialStarts = range(nPotentialStarts)
# 
#     if nReads.lower() == 'tile' or int(nReads) >= nPotentialStarts:
#         if OPTIONS.uniqueReadsOnly:
#             starts = potentialStarts
#         else:
#             starts = map( lambda x: random.choice(potentialStarts), range(nPotentialStarts))
#     else:
#         starts = random.sample(potentialStarts, nPotentialStarts)
# 
#     #print 'potentialStarts', potentialStarts
#     #print 'Starts', starts
#     
#     def sampleRead( start ):
#         if mutSeq[start] <> '-':
#             read, stop = accumulateBases( mutSeq, start, readLen )
#             remaining = len(mutSeq) - stop
#             refForRead = refSeq[start:stop]
#             #print 'XXXX', start, stop, refForRead, read
#             cigar = mutationsCIGAR( alignStart, readLen, mutations )
#             if not OPTIONS.quiet:
#                 leftPads = string.join(start * ['-'], '')
#                 rightPads = string.join(remaining * ['-'], '')
#                 print leftPads + mutSeq.tostring()[start:stop] + rightPads, start, stop, cigar
#             return filter(notDash, read), alignStart + start, cigar 
#         else:
#             return False, False, False
#     
#     reads = map( sampleRead, starts )
#     return [read for read in reads if read[0] <> False]

class refSite:
    def __init__( self, refSeq, mutSeq, alignStart, mutations ):
        self.refSeq = refSeq
        self.mutSeq = mutSeq
        self.alignStart = alignStart
        self.mutations = mutations

def sample1Read( refSite, start, readLen, reverseComplementP = False ):
    if refSite.mutSeq[start] <> '-':
        read, stop = accumulateBases( refSite.mutSeq, start, readLen )
        
        mutSubSeq = refSite.mutSeq[start:stop]
        if reverseComplementP:
            s = Seq( read, IUPAC.unambiguous_dna ) 
            #print 'reverseComplementP', read, s.reverse_complement().tostring(), mutSubSeq
            read = s.reverse_complement().tostring()
            mutSubSeq.complement()
        
        remaining = len(refSite.mutSeq) - stop
        refForRead = refSite.refSeq[start:stop]
        #print 'XXXX', start, stop, refForRead, read
        cigar = mutationsCIGAR( start, refSite.alignStart, readLen, refSite.mutations )
        leftPads = string.join(start * ['-'], '')
        rightPads = string.join(remaining * ['-'], '')
        ppReadStr = leftPads + mutSubSeq.tostring() + rightPads # + '(' + read + ')'
        return True, filter(notDash, read), refSite.alignStart + start, cigar, ppReadStr
    else:
        return False, None, None, None, None

def sampleReadsFromAlignment(wholeRef, refSeq, mutSeq, alignStart, readLen, nReads, mutations, pairedOffset, mutations2, refSeq2, mutSeq2):
    #print 'REF:', refSeq.tostring()
    #print 'MUT:', mutSeq.tostring()
    #print ''

    # pairedSite, mutations2, refSeq2, mutSeq2 can all be None
    
    refSite1 = refSite( refSeq, mutSeq, alignStart, mutations )
    refSite2 = refSite( refSeq2, mutSeq2, alignStart + pairedOffset, mutations2 )
    
    #print refSite1.alignStart, refSite2.alignStart, pairedOffset
    
    # we are assuming that mutSeq is exactly 1 + 2 * readLen in size, which the mutation
    # right in the middle
    lastGoodStartSite = readLen - 1
    if not OPTIONS.quiet:
        if OPTIONS.generatePairedEnds:
            r = wholeRef.seq[refSite1.alignStart:refSite2.alignStart + 2*readLen - 1]
            print r.tostring(), 'ref'
            print r.complement().tostring(), 'ref, complement'
        else:
            print refSeq.tostring(), 'ref'
            print mutSeq.tostring(), 'mut'
    
    # we can potentially start at any site in mutSeq
    nPotentialStarts = len(mutSeq) - readLen + 1
    potentialStarts = range(nPotentialStarts)
    #print 'potentialStarts', potentialStarts

    if nReads.lower() == 'tile' or int(nReads) >= nPotentialStarts:
        if OPTIONS.uniqueReadsOnly:
            starts = potentialStarts
        else:
            starts = map( lambda x: random.choice(potentialStarts), range(nPotentialStarts))
    else:
        starts = random.sample(potentialStarts, int(nReads))

    #print 'Starts', starts
   
    def sampleRead( start ):
        good1, read1, pos1, cigar1, ppstr1 = sample1Read( refSite1, start, readLen )
        
        if OPTIONS.generatePairedEnds:
            good2, read2, pos2, cigar2, ppstr2 = sample1Read( refSite2, start, readLen, True )
            if not OPTIONS.quiet:
                offsetDashes = string.join(['-'] * (pairedOffset), '') 
                print
                print ppstr1 + offsetDashes, pos1, cigar1
                print offsetDashes + ppstr2, pos2, cigar2
            return [(read1, pos1, cigar1), (read2, pos2, cigar2)]
        else:
            if not OPTIONS.quiet:
                print ppstr1, pos1, cigar1
            return [(read1, pos1, cigar1), (None, 0, None)]
    
    reads = map( sampleRead, starts )
    return [read for read in reads if read[0][0] <> None]

def notDash( c ):
    return c <> '-'

def fakeQuals( seq ):
    return len(seq) * [30]
    #return range(len(seq))

def alignedRead2SAM( siteID, readID, fastaID, read, pos, cigar, read2, pos2, cigar2 ):
    rname = fastaID
    mapq = 40
    quals = fakeQuals( read )

    qname = '%s:%d.read.%d' % (fastaID, siteID, readID)
    
    if not OPTIONS.generatePairedEnds:
        # not in paired mode
        flags = []
        record1 = SAMRecordFromArgs( qname, flags, rname, pos, mapq, cigar, read, quals )
        record2 = None
    else:
        flags = [SAM_SEQPAIRED, SAM_MAPPAIRED]
        insertSize = pos2 - pos - len(read)
        record1 = SAMRecordFromArgs( qname, flags, rname, pos, mapq, cigar, read, quals, pairContig = rname, pairPos = pos2, insertSize = insertSize )
        record2 = SAMRecordFromArgs( qname + 'p', flags, rname, pos2, mapq, cigar2, read2, quals, pairContig = rname, pairPos = pos, insertSize = insertSize )
        
        #print 'read1, read2', record1.getSeq(), record2.getSeq()
        
    return [record1, record2]

# ---------------------------------------------------------------------------
#
# paired end reads
#
# ---------------------------------------------------------------------------
def pairedReadSite( ref, leftMutSite, readLen ):
    pairedOffset = int(round(random.normalvariate(OPTIONS.insertSize, OPTIONS.insertDev)))

    def fail(): return False, None
    
    if pairedOffset < 0:
        return fail()
        
    pairedSite = pairedOffset + leftMutSite + readLen
    #print 'pairedSite', pairedOffset, leftMutSite, pairedSite
    
    if pairedSite + readLen >= len(ref):
        return fail()
    refsub = ref[pairedSite - readLen:pairedSite + readLen + 1]
    if isBadSequenceForSim( refsub ):
        return fail()
    
    return True, pairedSite
    


# ---------------------------------------------------------------------------
#
# build allowed regions and sampling starts from region
#
# ---------------------------------------------------------------------------
def parseSamplingRange( arg, ref, readLen ):
    # returns a list of the allowed regions to sample in the reference

    def one(x):
        elts = x.split()
        #print 'elts', elts
        if len(elts) > 0:
            elts1 = elts[0].split(':')
            #print 'elts1', elts1
            return map( int, elts1 )
        else:
            return False
  
    if arg <> None:
        try:
            return [ one(arg) ]
        except:
            # try to read it as a file
            return filter( None, map( one, open(arg).readlines() ) )
    else:
        return [[0, len(ref.seq) - readLen]]

def sampleStartSites( sampleRanges, nsites ):
    print 'sampleStartSites', sampleRanges, nsites
    
    # build a set of potential sites
    if len(sampleRanges) > 1:
        potentialSites = set()
        for start, end in sampleRanges:
            #print 'potentialSites', start, end, potentialSites
            potentialSites = potentialSites.union(set(xrange(start, end))) 
    else:
        start, end = sampleRanges[0]
        potentialSites = xrange(start, end) 

    #print 'potentialSites', potentialSites
    
    # choose sites from potentials
    if len(potentialSites) <= nsites:
        # tile every site
        sites = potentialSites
    else:
        # we need to sample from the potential sites
        sites = random.sample( potentialSites, nsites )

    # print 'Sites', sites
    print 'Number of potential start sites:', len(potentialSites)
    print 'Number of start sites that will be generated:', len(sites)
    return sites 

def readRef(referenceFasta):
	handle = open(referenceFasta)
	for seq_record in SeqIO.parse(handle, "fasta"):
		print seq_record.id 
		print seq_record.name
		#print repr(seq_record.seq)
		print len(seq_record.seq)
		yield seq_record
	handle.close()

import os.path
def outputFilename():
    root = os.path.splitext(os.path.split(OPTIONS.reference)[1])[0]
    if OPTIONS.generatePairedEnds:
        pairedP = 'Yes'
    else:
        pairedP = 'No'
        
    params = [['mut',OPTIONS.mutationType], ['den', max(OPTIONS.mutationDensity, OPTIONS.indelSize)], [OPTIONS.readLen,'bp'], \
              ['nsites', OPTIONS.nSites], ['cov', OPTIONS.coverage], ['range', OPTIONS.range], ['paired', pairedP]]
    filename = root + '__' + string.join(map( lambda p: string.join(map(str, p),'.'), params), '_')
    root = os.path.join(OPTIONS.outputPrefix, filename)
    return root, root + '.sam'

def main():
    global OPTIONS, ROOT

    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-r", "--ref", dest="reference",
                        type="string", default=None,
                        help="Reference DNA seqeunce in FASTA format")
    parser.add_option("-o", "--outputPrefix", dest="outputPrefix",
                        type="string", default='',
                        help="Output Prefix SAM file")
    parser.add_option("-c", "--coverage", dest="coverage",
                        type="string", default='1',
                        help="Number of reads to produce that cover each mutated region for the reference, or tile for all possible reads.")
    parser.add_option("-n", "--nsites", dest="nSites",
                        type=int, default=1,
                        help="Number of sites to generate mutations in the reference")
    parser.add_option("-l", "--readlen", dest="readLen",
                        type=int, default=36,
                        help="Length of reads to generate")
    parser.add_option("-s", "--seed", dest="seed",
                        type=int, default=123456789,
                        help="Random number seed.  If 0, then system clock is used")
    parser.add_option("-t", "--muttype", dest="mutationType",
                        type="string", default='None',
                        help="Type of mutation to introduce from reference.  Can be None, SNP, insertion, or deletion")
    parser.add_option("-e", "--mutdensity", dest="mutationDensity",
                        type=int, default=1,
                        help="Number of mutations to introduce from the reference")
    parser.add_option("-x", "--range", dest="range",
                        type="string", default=None,
                        help="Sampling range restriction, in the form of x:y without a space")
    parser.add_option("-u", "--unique", dest="uniqueReadsOnly",
                        action='store_true', default=True,
                        help="Assumes that the user wants unique reads generated.  If the coverage is greater than the number of unique reads at the site, the region is just tiled")
    parser.add_option("-i", "--indelsize", dest="indelSize",
                        type=int, default=0,
                        help="Size in bp of the insertion or deletion to introduce")

    # Paired end options
    parser.add_option("", "--paired", dest="generatePairedEnds",
                        action='store_true', default=False,
                        help="Should we generate paired end reads?")
    parser.add_option("", "--insertSize", dest="insertSize",
                        type=int, default=280,
                        help="Size in bp of the paired library insertion")
    parser.add_option("", "--insertDev", dest="insertDev",
                        type=int, default=5,
                        help="Standard deviation in the size in bp of the paired library insertion")

    parser.add_option("-q", "--quiet", dest="quiet",
                        action='store_true', default=False,
                        help="Be quiet when generating output")
    parser.add_option("-d", "--debug", dest="debug",
                        action='store_true', default=False,
                        help="Verbose debugging output?")
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 0:
        parser.error("incorrect number of arguments")
	
    root, outputSAM = outputFilename()

    if OPTIONS.seed == 0:
        random.seed(None)
    else:
        random.seed(OPTIONS.seed)

    OPTIONS.mutationType = OPTIONS.mutationType.upper()
    
    print '------------------------------------------------------------' 
    print 'program:', sys.argv[0]
    print ''
    print '  ref:               ', OPTIONS.reference
    print '  output:            ', outputSAM
    print ''
    print '  mutation type:     ', OPTIONS.mutationType
    print '  mutation density:  ', OPTIONS.mutationDensity
    print '  indel size:        ', OPTIONS.indelSize
    print '  readLen:           ', OPTIONS.readLen
    print '  nSites:            ', OPTIONS.nSites
    print '  coverage:          ', OPTIONS.coverage
    print '  range restriction: ', OPTIONS.range
    print ''
    print '  paired ends?:      ', OPTIONS.generatePairedEnds
    print '  insert size:       ', OPTIONS.insertSize
    print '  insert stddev:     ', OPTIONS.insertDev
    print ''
    print '  Debugging?:        ', OPTIONS.debug
    print '------------------------------------------------------------'

    if OPTIONS.mutationType <> 'SNP' and OPTIONS.mutationDensity > 1:
        raise Exception('Does not support mutation density > 1 for mutations of class', OPTIONS.mutationType)

    readLen = OPTIONS.readLen
    fastaRecords = [seq for seq in readRef(OPTIONS.reference)]
    header = SAMHeader( fastaRecords[0].id, len(fastaRecords[0].seq) ) 
    SAMout = SAMIO( outputSAM, header, debugging=OPTIONS.debug )
    
    mutationsout = open(root + '.mutations.txt', 'w')
    qualsout = open(root + '.quals.txt', 'w')
    refout = open(root + '.ref', 'w')

    fastaout = open(root + '.fasta', 'w')
    if OPTIONS.generatePairedEnds:
        fastaout2 = open(root + '.pairs.fasta', 'w')
        qualsout2 = open(root + '.pairs.quals.txt', 'w')

    counter = 0
    refLen = 0
    for ref in fastaRecords:
        refLen = len(ref.seq)
        
        # write the crazy ref file info needed by samtools
        print >> refout, ref.id + '\t' + str(refLen)

        sampleRanges = parseSamplingRange( OPTIONS.range, ref, readLen )
        sites = sampleStartSites( sampleRanges, OPTIONS.nSites )
        
        for mutSite, siteCounter in zip( sites, range(1, len(sites)+1) ):
            #print 'Sampling site:', mutSite, refLen 
            #mutSite = readLen-1
            nSitesSelected = max(int(round(len(sites)/100.0)), 1)
            #print siteCounter, len(sites), nSitesSelected)
            if siteCounter % nSitesSelected == 0:
                print 'Sampled site %d of %d (%.2f%%)' % ( siteCounter, len(sites), (100.0*siteCounter) / len(sites)) 
                            
            good, mutations, refSeq, mutSeq = mutateReference(ref.seq, ref.id, mutSite, OPTIONS.mutationType, [], readLen)

            pairedSite, mutations2, refSeq2, mutSeq2 = 0, None, None, None
            if good and OPTIONS.generatePairedEnds:
                if OPTIONS.mutationType == 'SNP':
                    pairedMutType = 'SNP'
                else:
                    pairedMutType = 'NONE'
                
                good, pairedSite = pairedReadSite( ref.seq, mutSite, readLen )
                print 'pairedReadSite', good, pairedSite, good
                if good:
                    good, mutations2, refSeq2, mutSeq2 = mutateReference(ref.seq, ref.id, pairedSite, pairedMutType, [], readLen)

            #print 'Good', good, mutations, refSeq, mutSeq
            if good:
                print >> mutationsout, ref.id + ':' + str(mutSite) + '\t' + string.join(map(str, mutations), '\t')
                #print 'Mutations', mutations

                refSeqPos = mutSite - readLen + 1
                readPairs = sampleReadsFromAlignment(ref, refSeq, mutSeq, refSeqPos, readLen, OPTIONS.coverage, mutations, pairedSite - mutSite, mutations2, refSeq2, mutSeq2)
                for i in range(len(readPairs)):
                    counter += 1

                    read1, read2 = readPairs[i]
                    seq, pos, cigar = read1
                    seq2, pos2, cigar2 = read2

                    # write out the sam and fasta files
                    def write1( record, recordi ):
                        if record <> None:
                            SAMout.writeRecord( record )
                            sequences = [SeqRecord(Seq(record.getSeq()), id = record.getQName(), description='')]
                            #print sequences
                            
                            if recordi % 2 == 0:
                                localFastaOut, localQualsOut = fastaout, qualsout
                            else:
                                localFastaOut, localQualsOut = fastaout2, qualsout2

                            SeqIO.write(sequences, localFastaOut, "fasta")
                            print >> localQualsOut, record.getQName(), string.join(map( str, record.getQuals() ), ' ')
                    
                    #print i, read1, read2
                    records = alignedRead2SAM( mutSite, i+1, ref.id, seq, pos + 1, cigar, seq2, pos2 + 1, cigar2 )
                    map( write1, records, range( len(records) ) )

    SAMout.close()
    qualsout.close()
    fastaout.close()
    refout.close()
    mutationsout.close()
    
    if OPTIONS.generatePairedEnds:
        qualsout2.close()
        fastaout2.close()
    
#    for record in SAMIO( outputSAM ):
#        print record
	
if __name__ == "__main__":
    import Bio
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC, Gapped
    main()


