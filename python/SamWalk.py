#!/util/bin/python
#
# 
# This file really needs to be cleaned up and functions need to be improved
#
# (1) utility operations should be moved to a utility file
# (2) Need to add support for windows paired reads, information about where
#     reads come from (which file), and n-base pair loci
# (3) Move to C++ based SAM/BAM i/o system
# (4) Figure out a way to get the reference to stream, and use a automatically pickled version
# (5) Add capability to walk over all reference bases, regardless of coverage.  Currently
#     we only walk over covered bases.
#
#
from optparse import OptionParser
import os
import sys
import random
import string
import os.path
import heapq
from itertools import *

from SAM import *
import SimpleSAM
#import SimulateReads    # utility functions should be imported

# Global variable containing command line arguments
OPTIONS = None

#
# Trivial class that added pushback support to a stream
#
class peakStream(object):
    """A simple wrapper around a stream that adds unget operation"""
    def __init__( self, stream ):
        self.ungot = []
        self.s = stream
        
    def unget( self, elt ):
        """Put elt at the front of the underlying stream, making the next get call return it"""
        
        #print 'unget', self.ungot
        self.ungot.append(elt)
        #print '  -> unget', self.ungot

    # we hold the next() operation
    def __iter__( self ): return self

    def next(self):
        """For python iterator support"""
        #print 'ungot', self.ungot
        if self.ungot <> []:
            elt = self.ungot.pop(0)
        else:
            elt = self.s.next()

        #print '  -> ungot', self.ungot
        #print '  elt', elt
        return elt

# 
# Returns a single, combined stream of SAM records in chromosome order from a list of 
# input streams, each in reference order
# 
def rawReadStream( ref, inputs ):
    """Returns a single iterable that returns records in reference order
    from each of the input streams"""
    def includePriorities( stream ):
        return imap( lambda x: (x.getPos(), x), stream )
    def prunePriorities( stream ):
        return imap( lambda p: SimpleSAM.MappedReadFromSamRecord(ref, p[1]), stream )

    with_priorities = map( includePriorities, inputs )
    return prunePriorities( heapq.merge( *with_priorities ) )

# 
# Just wraps the raw read stream objects in a list as they come out
# 
def readStream( ref, inputs ):
    """Returns a single iterable that returns SAM records in reference order
    from each of the input streams"""
    rs = rawReadStream( ref, inputs )
    return imap( lambda x: [x], rs )

#
# More complex stream object.  Takes a list of input streams and creates a stream
# returning successive loci covered by the reads in the combined input stream.
#
class lociStream():
    """A streaming data structure that returns reads spanning each covered loci in
    the input reads, offsets into them
    where the bases are equivalent, and the position of the locus.
    """
    def __init__( self, ref, inputs ):
        self.ref = ref
        self.rs = peakStream(rawReadStream( ref, inputs ))
        self.window = []
        
        self.chr = None
        self.pos = None
        
    def __iter__(self):
        return self
        
    def pushRead( self, read ):
        self.window.append(read)
        
    def cleanWindow( self, chr, pos ):
        """Walk over the window of reads, deleting those who no longer cover the current pos"""
        if OPTIONS.debug:
            print 'cleanWindow start:', chr, pos, self.window 

        def keepReadP( read ):
            return read.chr == chr and pos >= read.start and pos <= read.end
        self.window = filter( keepReadP, self.window )

        if OPTIONS.debug:
            print 'cleanWindow stop:', chr, pos, self.window 

    def expandWindow( self, chr, pos ):
        """Keep reading from the read stream until we've added all reads from the stream covering pos"""
        if OPTIONS.debug:
            print 'expandWindow start:', chr, pos, self.window 
        for read in self.rs:
            #print 'read', read, pos
            if read.chr == chr and read.start <= pos and read.end >= pos:
                self.pushRead(read)
            else:
                self.rs.unget( read )
                #self.rs = chain( [read], self.rs )
                break
        if OPTIONS.debug:
            print 'expandWindow stop:', chr, pos, self.window 
            
            
    #
    # This is the workhorse routine.  It constructs a window of reads covering the
    # next locus in the genome, and returns the reference, the contig, its position
    # in the genome, along with a vector of offsets into the list of reads that denote
    # the equivalent base among all the reads and reference.
    #
    def next(self):
        if self.pos <> None:
            self.cleanWindow(self.chr, self.pos)             # consume reads not covering pos
            self.expandWindow(self.chr, self.pos)            # add reads to window covering pos

        if self.window == []:
            # the window is empty, we need to jump to the first pos of the first read in the stream:
            nextRead = self.rs.next()
            self.pushRead( nextRead )
            self.chr = nextRead.chr
            self.pos = nextRead.start
            return self.next()
        else:
            # at this point, window contains all reads covering the pos, we need to return them 
            # and the offsets into each read for this loci
            def calcOffset( read ):
                offset = self.pos - read.start
                return offset
            
            offsets = map(calcOffset, self.window)
            currPos = self.pos
            self.pos += 1       # we are going to try to get the next position
            return self.ref, self.chr, currPos, offsets, self.window
            
# 
# Reference reader
#
def readRef(referenceFasta):
    ref = {}

    header = None
    seq = ''
    for line in open(referenceFasta):
        if line[0] == '>':
            if header <> None:
                ref[header] = seq.lower()
            seq = ''
            header = line[1:].strip()
        else:
            seq += line.strip()

    ref[header] = seq.lower()
    #print ref
    return ref

# 
# Main() procedure
#
def main():
    global OPTIONS, ROOT

    # ------------------------------------------------------------
    #
    # Setup command line options
    #
    # ------------------------------------------------------------
    usage = "usage: %prog [options] Walker [sam or bam file or file list]+"
    parser = OptionParser(usage=usage)
    parser.add_option("-r", "--ref", dest="reference",
                        type="string", default=None,
                        help="Reference DNA seqeunce in FASTA format")
    parser.add_option("-m", "--maxRecords", dest="maxRecords",
                        type=int, default=None,
                        help="Max number of SAM records to process")
    parser.add_option("-q", "--quiet", dest="quiet",
                        action='store_true', default=False,
                        help="Be quiet when generating output")
    parser.add_option("-d", "--debug", dest="debug",
                        action='store_true', default=False,
                        help="Verbose debugging output?")

    (OPTIONS, args) = parser.parse_args()
    print 'args', args
    if len(args) < 2:
        parser.error("Incorrect number of arguments")
	    

    # ------------------------------------------------------------
    #
    # Dynamically load the walker class (the first cmd line arg)
    # and initialize a walker class object
    # 
    # 
    # load walkers from standard location 
    #
    # ------------------------------------------------------------
    walkerName = args[0]
    walkerModule = __import__(walkerName, globals(), locals(), [], -1)
    walkerClass = walkerModule.__dict__[walkerName]
    walker = walkerClass()
    
    print walkerName
    print walkerModule
    print walkerClass
    print walker

    print '------------------------------------------------------------' 
    print 'program:', sys.argv[0]
    print ''
    print '  ref:               ', OPTIONS.reference
    print '  walker:            ', walker
    if walker.hasOption('name'):
        print '    ->               ', walker.getOption('name')
    if walker.hasOption('desc'):
        print '    ->               ', walker.getOption('desc')
    print '------------------------------------------------------------'


    # ------------------------------------------------------------
    #
    # Initialize the engine
    #
    # ------------------------------------------------------------
    
    # read the reference
    refs = readRef(OPTIONS.reference)

    # create the low-level SAMIO streams
    files = args[1:]
    inputs = map( lambda file: SAMIO( file, debugging=OPTIONS.debug ), files )
    
    # build the higher level stream object, either by record or by loci
    if walker.isRecordWalker():
        stream = readStream( refs, inputs )
    if walker.isLociWalker():
        stream = lociStream( refs, inputs )

    # ------------------------------------------------------------
    #
    # Move the walker object over the stream
    # 
    # For each element in the record or loci stream, invoke the walker 
    # object on it.  Use filter, map, and reduce to construct the sum
    #
    # ------------------------------------------------------------
    sum = walker.reduceDefault
    counter = 0
    for elt in stream:
        counter += 1
        #print 'processing elt', counter, elt, OPTIONS.maxRecords
        if OPTIONS.maxRecords <> None:
            if (counter * 10) % OPTIONS.maxRecords == 0:
                print '  => ', (counter * 100.0) / OPTIONS.maxRecords, '% done' 
        
            if counter > OPTIONS.maxRecords:
                break
        if walker.filterfunc( *elt ):
            x = walker.mapfunc( *elt )
            sum = walker.reducefunc( x, sum )
    print 'Result is', sum
    
if __name__ == "__main__":
    if True:
        import cProfile, pstats
        cProfile.run('main()', 'prof')
        stats = pstats.Stats('prof')
        stats.sort_stats('time')
        stats.print_stats(20)
    else:
        main()


