#!/util/bin/python

import Walker
from SimpleSAM import *

class WalkLociTest(object, Walker.Walker):
    def __init__(self):
        Walker.Walker.__init__( self, 'byLoci', {}, name = 'WalkLociTest' )
        self.debug = False

    def filterfunc( self, *data ):
        return True
        #return datum.getPos() % 2 == 0

    def mapfunc( self, ref, chrom, pos, offsets, reads ):

        if self.debug:
            print '------------------------------'
            print '  locus', chrom, pos
            print '  offsets = ', offsets
            print '  reads = ', reads
            
        # 
        # This is the heart of the test walker.  You've got the reference, contig, 
        # pos, offsets, and reads.  The line below generate pile-ups
        #
        def getField(field):
            return map( lambda r, o: r.__getattribute__(field)[o], reads, offsets ) 

        def toASCII33( quals ):
            return ''.join( map( lambda q: chr(q+33), quals ))

        if False:
            # actually do some useful work
            refbase = ref[chrom][pos-1]
            bases = ''.join(getField('bases'))
            quals = toASCII33(getField('quals'))
            print chrom, pos, refbase, len(reads), all(map(lambda b: b == refbase, bases)), bases, quals
            
        if self.debug:
            for offset, read in zip( offsets, reads ):
                if self.debug:
                    print '    offset, read', offset, read
                print read.bases[offset], read.quals[offset]
        
        #print datum
        return 1
        #return len(datum.getSeq())
        #return datum.getSeq() + '/'
        #return "\n>%s\n%s" % (datum.getQName(), datum.getSeq())
     
    #reduceDefault = ''
    def reducefunc( self, x, sum ):
        #print x
        return x + sum



