package org.broadinstitute.sting.playground.gatk.walkers.graphalign;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.fasta.FastaReferenceWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.graph.SimpleDirectedGraph;
import org.apache.log4j.Logger;

import java.util.*;
import java.io.*;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.util.StringUtil;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

/**
 * A completely experimental read walker that consumes a graphical reference emitted by GraphReferenceBuilder as a
 * serialized java object and evaluates the number of mismatches to both the flat reference and the graphical
 * reference for each read [Not for public use and will change drastically in the future].
 */
public class GraphReferenceAssessor extends ReadWalker<Integer, Integer> {
    @Argument(fullName="graphFile", shortName="GF", doc="", required=true)
    String graphFile = null;
    ObjectInputStream graphSerialStream = null;

    @Argument(fullName="MAX", shortName="MAX", doc="", required=false)
    int MAXREADS = -1;

    @Argument(fullName="ignore0MM", shortName="I0", doc="", required=false)
    boolean IGNORE_0_MM = false;

    @Argument(fullName="DEBUG", shortName="DB", doc="", required=false)
    int DEBUG_LEVEL = 0;

    static boolean DEBUG = false;
    static boolean DEBUG2 = false; // crazy level

    @Argument(fullName="read", doc="", required=false)
    String onlyDoRead = null;

    ReferenceGraph graphRef = null;

    public void initialize() {
        super.initialize();

        DEBUG = DEBUG_LEVEL > 0;
        DEBUG2 = DEBUG_LEVEL > 1; // crazy level
        
        try {
            logger.info("Reading graph reference " + graphFile );
            graphSerialStream = new ObjectInputStream( new FileInputStream( graphFile ) );
            graphRef = (ReferenceGraph)graphSerialStream.readObject();
            graphRef.setDebugPrinting(DEBUG);
            graphRef.validateGraph();            
            logger.info(graphRef.toBriefString());
        } catch ( FileNotFoundException e ) {
            throw new StingException("Couldn't open file " + graphFile, e);
        } catch ( IOException e ) {
            throw new StingException("Couldn't write to file " + graphFile, e);
        } catch ( ClassNotFoundException e ) {
            throw new StingException("Couldn't read ReferenceGraph from file " + graphFile, e);
        }
    }

    private static MismatchCounter countMismatches(byte[] ref, int refOffset, byte[] bases, byte[] quals, int basesOffset, int length) {
        MismatchCounter mm = new MismatchCounter();

        for ( int i = 0; i < length; i++ ) {
            byte rawRefBase = ref[i + refOffset];
            byte rawReadBase = bases[i + basesOffset];
            int fragBase = BaseUtils.simpleBaseToBaseIndex((char)rawRefBase);
            int readBase = BaseUtils.simpleBaseToBaseIndex((char)rawReadBase);

            boolean mmP = fragBase != -1 && readBase != -1 && fragBase != readBase;
            if ( mmP ) {
                mm.nMM++;
                mm.qSum += quals != null ? quals[i + basesOffset] : 0;
            }

            if ( GraphReferenceAssessor.DEBUG2 )
                System.out.printf("%s%d %c %c %s %b%n", Utils.dupString(' ', basesOffset + 2), basesOffset, (char)rawRefBase, (char)rawReadBase, mm, mmP);
        }

        return mm;
    }

    private static MismatchCounter countMismatches(byte[] ref, byte[] bases, byte[] quals) {
        return countMismatches(ref, 0, bases, quals, 0, bases.length);
    }

    private static MismatchCounter countMismatchesOnGraph( ReferenceGraph graph, Collection<Fragment> frags, int fragOffset, byte[] bases, byte[] quals, int readOffset ) {
        if ( frags.size() == 0 )
            throw new RuntimeException("Fragment list is empty!");

        MismatchCounter minNMM = MismatchCounter.MAX_VALUE;

        for ( Fragment next : frags ) {
            MismatchCounter recNMM = countMismatchesOnGraph( graph, next, 0, bases, quals, readOffset );
            minNMM = minNMM.min( recNMM );
        }

        return minNMM;
    }

    private static MismatchCounter countMismatchesOnGraph( ReferenceGraph graph, Fragment frag, int fragOffset, byte[] bases, byte[] quals, int readOffset ) {
        if ( GraphReferenceAssessor.DEBUG )System.out.printf("%sfrag %s -> %d%n", Utils.dupString(' ', readOffset + 2), frag, readOffset);

        MismatchCounter mm = new MismatchCounter();

        if ( readOffset < bases.length ) {
            int nRemainingBases = bases.length - readOffset;
            int cmpLength = frag.getBaseLengthFrom(fragOffset, nRemainingBases);       // how many bases over in the fragment are we from the offset
            MismatchCounter fragMM = countMismatches(frag.getUnderlyingBases(), frag.getUnderlyingOffset() + fragOffset, bases, quals, readOffset, cmpLength);
            mm.add(fragMM);

//            // still have some counting to do
//            for ( int i = 0; i < baseLength; i++ ) {
//                int fragBaseOffset = fragOffset + i;
//                int readBaseOffset = readOffset + i;
//
//                byte rawFragBase = frag.getBase(fragBaseOffset);
//                byte rawReadBase = bases[readBaseOffset];
//                int fragBase = BaseUtils.simpleBaseToBaseIndex((char)rawFragBase);
//                int readBase = BaseUtils.simpleBaseToBaseIndex((char)rawReadBase);
//
//                boolean mmP = fragBase != -1 && readBase != -1 && fragBase != readBase;
//                if ( mmP ) nMM++;

            if ( nRemainingBases > cmpLength ) {
                MismatchCounter recMM = countMismatchesOnGraph( graph, graph.outgoingFragments(frag), 0, bases, quals, readOffset + cmpLength );
                mm.add(recMM);
            }
        }

        if ( GraphReferenceAssessor.DEBUG ) System.out.printf("%s=> %s%n", Utils.dupString(' ', readOffset + 2), mm);
        return mm;
    }

    private static MismatchCounter countMismatchesOnGraph(ReferenceGraph graph, SAMRecord read) {
        if ( GraphReferenceAssessor.DEBUG ) System.out.printf("countMismatchesOnGraph( read=%s%n", read.getReadName());
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(read);
        MismatchCounter minNMM = MismatchCounter.MAX_VALUE;

        for ( Fragment frag : graph.getStartingFragment(loc) ) {
            int fragOffset = frag.getFragOffsetFrom(loc);       // how many bases over in the fragment are we from the offset

            if ( GraphReferenceAssessor.DEBUG )
                System.out.printf("  countMismatchesOnGraph frag=%s loc=%s bases=%s offset=%d%n", frag, loc, read.getReadString(), fragOffset);

            MismatchCounter recNMM = countMismatchesOnGraph(graph, frag, fragOffset, read.getReadBases(), read.getBaseQualities(), 0);
            minNMM = minNMM.min( recNMM );
        }

        return minNMM;
    }

    public Integer map(char[] refArg, SAMRecord read, ReadMetaDataTracker metaDataTracker) {

        if ( MAXREADS-- == 0 ) {
            System.exit(0);
        } else if ( onlyDoRead != null && ! read.getReadName().equals(onlyDoRead) ) {
            ;
        } else if ( ! read.getReadUnmappedFlag() && read.getCigar().numCigarElements() == 1 ) {
            try {
                byte[] ref = BaseUtils.charSeq2byteSeq(refArg);
                // we're all XM
                int nMMFromRead = (Short)read.getAttribute("NM");
                MismatchCounter nAlignedMM = countMismatches(ref, read.getReadBases(), read.getBaseQualities());
                if ( ! IGNORE_0_MM || nAlignedMM.nMM > 0 ) {
                    MismatchCounter nGraphMM = countMismatchesOnGraph(graphRef, read);
                    MismatchCounter deltaMM = nAlignedMM.minus(nGraphMM);

                    out.printf("%50s with %5s at %10s: mismatches: %3d (delta %3d) -- %3d %3d -- %3d %3d -- delta %3d %3d%n",
                            read.getReadName(), read.getCigarString(), GenomeLocParser.createGenomeLoc(read),
                            nMMFromRead, nMMFromRead - nAlignedMM.nMM,
                            nAlignedMM.nMM, nAlignedMM.qSum,
                            nGraphMM.nMM, nGraphMM.qSum,
                            deltaMM.nMM, deltaMM.qSum);

                    if ( deltaMM.nMM < 0 || deltaMM.qSum < 0 )
                        throw new StingException(read.getReadName() + " is miscalculated");
                }
            } catch ( Exception e ) {
                System.out.printf("Exception at %s at %s%n", read.getReadName(), GenomeLocParser.createGenomeLoc(read));
                throw new RuntimeException(e);
            }
        } else {
                ; // don't do anything
        }

        return 0;
    }

    /**
     * reduceInit is called once before any calls to the map function.  We use it here to setup the output
     * bam file, if it was specified on the command line
     *
     * @return
     */
    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer readScore, Integer data) {
        return data + readScore;
    }

    public void onTraversalDone(Integer data) {
        out.printf(data.toString());
    }
}

class MismatchCounter {
    int nMM = 0;
    int qSum = 0;

    public static MismatchCounter MAX_VALUE = new MismatchCounter(Integer.MAX_VALUE, Integer.MAX_VALUE );

    public MismatchCounter() {}

    public MismatchCounter(int nMM, int qSum) {
        this.nMM = nMM;
        this.qSum = qSum;
    }

    public void add(MismatchCounter that) {
        this.nMM += that.nMM;
        this.qSum += that.qSum;
    }

    public MismatchCounter min(MismatchCounter that) {
        int cmpQSum = Integer.valueOf(this.qSum).compareTo(that.qSum);
        if ( cmpQSum < 0 ) { return this; }
        else if ( cmpQSum > 0 ) { return that; }
        else if ( this.nMM < that.nMM ) { return this; }
        else if ( this.nMM > that.nMM ) { return that; }
        else { return this; }
    }

    public MismatchCounter minus(MismatchCounter that) {
        return new MismatchCounter(this.nMM - that.nMM, this.qSum - that.qSum);
    }

    public String toString() { return String.format("[MM %d %d]", nMM, qSum); }
}