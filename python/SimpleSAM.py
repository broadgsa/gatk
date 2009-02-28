from SAM import *

def read2align( refSeq, readSeq, chr, pos, cigar ):
    """Converts the read bases, a reference sequence, and a cigar string into a pair of
    alignments strings for refSeq and readSeq.  The alignment string has base equivalents
    between the reference and the read bases"""
    # 
    # only works for non-indel containing reads right now!
    #
    contig = refSeq[chr]
    refAlign = contig[pos:pos+len(readSeq)]
    readAlign = readSeq
    return refAlign, readAlign

def MappedReadFromSamRecord(ref, r):
    """Converts SameRead into a SimpleSamRecord object"""
    #
    # This is a temp. function while we wait for Ted's C++ SAM/BAM io system to come online    
    #
    refAlign, readAlign = read2align(ref, r.getSeq(), r.getRname(), r.getPos(), r.getCigar)
    return MappedRead( r.getQName(), r.getFlags(), r.getRname(), r.getPos(), r.getMapq(), 
                       r.getCigar(), r.getSeq(), refAlign, readAlign, r.getQuals() ) 

class MappedRead(object):
    """Higher level object representing a mapped read.  
    
    Attempts to be more accessible than a raw SAM record, doing a lot of processing on the record w.r.t.
    the reference to make analysis of the reads easier.  This is the fundamental data structure
    used to represent reads through the walker engine"""
    
    def __init__( self, qname, flags, chr, pos, mapq, cigar, readBases, refAlign, readAlign, qualsByBase ):
        # basic information about the read
        self.qname = qname
        self.mapq = mapq
        self.cigar = cigar
        self.bases = readBases          # actual bases of reads, unaligned, may be fw or reverse
        self.quals = qualsByBase        # as a vector of floats, by base position

        # dealing with position information, more than a standard record Hasbro
        self.chr     = chr
        self.start   = pos
        self.readlen = len(readBases)
        self.end     = pos + len(refAlign) - 1

        # print 'refAlign', refAlign, len(refAlign)

        # dealing with stand
        self.flags = flags
        self.isForwardStrand = SAMFlagIsSet(flags, SAM_QUERYSTRAND)
        self.isReverseStrand = not SAMFlagIsSet(flags, SAM_QUERYSTRAND)

        self.isPrimary   = SAMFlagIsSet(flags, SAM_NOTPRIMARY)
        self.isUnmapped  = SAMFlagIsSet(flags, SAM_UNMAPPED)
        self.isPaired    = SAMFlagIsSet(flags, SAM_SEQPAIRED)
        self.isMapPaired = SAMFlagIsSet(flags, SAM_MAPPAIRED)
        self.isFirstReadInPair  = SAMFlagIsSet(flags, SAM_ISFIRSTREAD)
        self.isSecondReadInPair = SAMFlagIsSet(flags, SAM_ISSECONDREAD)

        #self.isMapPaired = SAMFlagIsSet(flags, SAM_MATEUNMAPPED)
        #self.isMapPaired = SAMFlagIsSet(flags, SAM_MATESTRAND)

        # everything will be reverse complemented already
        self.refAlign = refAlign        # 
        self.readAlign = readAlign      # all should be forward strand

    def __str__(self):
        return '<MappedRead %s:%d-%d>' % (self.chr, self.start, self.end)

    __repr__ = __str__ 