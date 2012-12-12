/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.iterators;

import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.downsampling.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.*;

/**
 * Iterator that traverses a SAM File, accumulating information on a per-locus basis
 */
public class LocusIteratorByState extends LocusIterator {
    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(LegacyLocusIteratorByState.class);

    // -----------------------------------------------------------------------------------------------------------------
    //
    // member fields
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Used to create new GenomeLocs.
     */
    private final GenomeLocParser genomeLocParser;
    private final ArrayList<String> samples;
    private final ReadStateManager readStates;

    protected static class SAMRecordState {
        SAMRecord read;
        int readOffset = -1;                    // how far are we offset from the start of the read bases?
        int genomeOffset = -1;                  // how far are we offset from the alignment start on the genome?

        Cigar cigar = null;
        int cigarOffset = -1;
        CigarElement curElement = null;
        int nCigarElements = 0;

        int cigarElementCounter = -1;           // how far are we into a single cigarElement

        // The logical model for generating extended events is as follows: the "record state" implements the traversal
        // along the reference; thus stepForwardOnGenome() returns on every and only on actual reference bases. This
        // can be a (mis)match or a deletion (in the latter case, we still return on every individual reference base the
        // deletion spans). In the extended events mode, the record state also remembers if there was an insertion, or
        // if the deletion just started *right before* the current reference base the record state is pointing to upon the return from
        // stepForwardOnGenome(). The next call to stepForwardOnGenome() will clear that memory (as we remember only extended
        // events immediately preceding the current reference base).

        public SAMRecordState(SAMRecord read) {
            this.read = read;
            cigar = read.getCigar();
            nCigarElements = cigar.numCigarElements();

            //System.out.printf("Creating a SAMRecordState: %s%n", this);
        }

        public SAMRecord getRead() {
            return read;
        }

        /**
         * What is our current offset in the read's bases that aligns us with the reference genome?
         *
         * @return
         */
        public int getReadOffset() {
            return readOffset;
        }

        /**
         * What is the current offset w.r.t. the alignment state that aligns us to the readOffset?
         *
         * @return
         */
        public int getGenomeOffset() {
            return genomeOffset;
        }

        public int getGenomePosition() {
            return read.getAlignmentStart() + getGenomeOffset();
        }

        public GenomeLoc getLocation(GenomeLocParser genomeLocParser) {
            return genomeLocParser.createGenomeLoc(read.getReferenceName(), getGenomePosition());
        }

        public CigarOperator getCurrentCigarOperator() {
            return curElement.getOperator();
        }

        public String toString() {
            return String.format("%s ro=%d go=%d co=%d cec=%d %s", read.getReadName(), readOffset, genomeOffset, cigarOffset, cigarElementCounter, curElement);
        }

        public CigarElement peekForwardOnGenome() {
            return ( cigarElementCounter + 1 > curElement.getLength() && cigarOffset + 1 < nCigarElements ? cigar.getCigarElement(cigarOffset + 1) : curElement );
        }

        public CigarElement peekBackwardOnGenome() {
            return ( cigarElementCounter - 1 == 0 && cigarOffset - 1 > 0 ? cigar.getCigarElement(cigarOffset - 1) : curElement );
        }

        
        public CigarOperator stepForwardOnGenome() {
            // we enter this method with readOffset = index of the last processed base on the read
            // (-1 if we did not process a single base yet); this can be last matching base, or last base of an insertion


            if (curElement == null || ++cigarElementCounter > curElement.getLength()) {
                cigarOffset++;
                if (cigarOffset < nCigarElements) {
                    curElement = cigar.getCigarElement(cigarOffset);
                    cigarElementCounter = 0;
                    // next line: guards against cigar elements of length 0; when new cigar element is retrieved,
                    // we reenter in order to re-check cigarElementCounter against curElement's length
                    return stepForwardOnGenome();
                } else {
                    if (curElement != null && curElement.getOperator() == CigarOperator.D)
                        throw new UserException.MalformedBAM(read, "read ends with deletion. Cigar: " + read.getCigarString() + ". Although the SAM spec technically permits such reads, this is often indicative of malformed files. If you are sure you want to use this file, re-run your analysis with the extra option: -rf BadCigar");
                        
                    // Reads that contain indels model the genomeOffset as the following base in the reference.  Because
                    // we fall into this else block only when indels end the read, increment genomeOffset  such that the
                    // current offset of this read is the next ref base after the end of the indel.  This position will
                    // model a point on the reference somewhere after the end of the read.
                    genomeOffset++; // extended events need that. Logically, it's legal to advance the genomic offset here:
                    // we do step forward on the ref, and by returning null we also indicate that we are past the read end.

                    return null;
                }
            }

            boolean done = false;
            switch (curElement.getOperator()) {
                case H: // ignore hard clips
                case P: // ignore pads
                    cigarElementCounter = curElement.getLength();
                    break;
                case I: // insertion w.r.t. the reference
                case S: // soft clip
                    cigarElementCounter = curElement.getLength();
                    readOffset += curElement.getLength();
                    break;
                case D: // deletion w.r.t. the reference
                    if (readOffset < 0)             // we don't want reads starting with deletion, this is a malformed cigar string
                        throw new UserException.MalformedBAM(read, "read starts with deletion. Cigar: " + read.getCigarString() + ". Although the SAM spec technically permits such reads, this is often indicative of malformed files. If you are sure you want to use this file, re-run your analysis with the extra option: -rf BadCigar");
                    // should be the same as N case
                    genomeOffset++;
                    done = true;
                    break;
                case N: // reference skip (looks and gets processed just like a "deletion", just different logical meaning)
                    genomeOffset++;
                    done = true;
                    break;
                case M:
                case EQ:
                case X:
                    readOffset++;
                    genomeOffset++;
                    done = true;
                    break;
                default:
                    throw new IllegalStateException("Case statement didn't deal with cigar op: " + curElement.getOperator());
            }

            return done ? curElement.getOperator() : stepForwardOnGenome();
        }
    }

    //final boolean DEBUG = false;
    //final boolean DEBUG2 = false && DEBUG;
    private ReadProperties readInfo;
    private AlignmentContext nextAlignmentContext;
    private boolean performDownsampling;

    // -----------------------------------------------------------------------------------------------------------------
    //
    // constructors and other basic operations
    //
    // -----------------------------------------------------------------------------------------------------------------

    public LocusIteratorByState(final Iterator<SAMRecord> samIterator, ReadProperties readInformation, GenomeLocParser genomeLocParser, Collection<String> samples) {
        this.readInfo = readInformation;
        this.genomeLocParser = genomeLocParser;
        this.samples = new ArrayList<String>(samples);

        // LIBS will invoke the Reservoir and Leveling downsamplers on the read stream if we're
        // downsampling to coverage by sample. SAMDataSource will have refrained from applying
        // any downsamplers to the read stream in this case, in the expectation that LIBS will
        // manage the downsampling. The reason for this is twofold: performance (don't have to
        // split/re-assemble the read stream in SAMDataSource), and to enable partial downsampling
        // of reads (eg., using half of a read, and throwing the rest away).
        this.performDownsampling = readInfo.getDownsamplingMethod() != null &&
                                   readInfo.getDownsamplingMethod().type == DownsampleType.BY_SAMPLE &&
                                   readInfo.getDownsamplingMethod().toCoverage != null;

        this.readStates = new ReadStateManager(samIterator);

        // currently the GATK expects this LocusIteratorByState to accept empty sample lists, when
        // there's no read data.  So we need to throw this error only when samIterator.hasNext() is true
        if (this.samples.isEmpty() && samIterator.hasNext()) {
            throw new IllegalArgumentException("samples list must not be empty");
        }
    }

    /**
     * For testing only.  Assumes that the incoming SAMRecords have no read groups, so creates a dummy sample list
     * for the system.
     */
    public final static Collection<String> sampleListForSAMWithoutReadGroups() {
        List<String> samples = new ArrayList<String>();
        samples.add(null);
        return samples;
    }

    public Iterator<AlignmentContext> iterator() {
        return this;
    }

    public void close() {
        //this.it.close();
    }

    public boolean hasNext() {
        lazyLoadNextAlignmentContext();
        return (nextAlignmentContext != null);
        //if ( DEBUG ) System.out.printf("hasNext() = %b%n", r);
    }

    private GenomeLoc getLocation() {
        return readStates.isEmpty() ? null : readStates.getFirst().getLocation(genomeLocParser);
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    // next() routine and associated collection operations
    //
    // -----------------------------------------------------------------------------------------------------------------
    public AlignmentContext next() {
        lazyLoadNextAlignmentContext();
        if (!hasNext())
            throw new NoSuchElementException("LocusIteratorByState: out of elements.");
        AlignmentContext currentAlignmentContext = nextAlignmentContext;
        nextAlignmentContext = null;
        return currentAlignmentContext;
    }

    /**
     * Creates the next alignment context from the given state.  Note that this is implemented as a lazy load method.
     * nextAlignmentContext MUST BE null in order for this method to advance to the next entry.
     */
    private void lazyLoadNextAlignmentContext() {
        while (nextAlignmentContext == null && readStates.hasNext()) {
            readStates.collectPendingReads();

            final GenomeLoc location = getLocation();
            final Map<String, ReadBackedPileupImpl> fullPileup = new HashMap<String, ReadBackedPileupImpl>();

            // TODO: How can you determine here whether the current pileup has been downsampled?
            boolean hasBeenSampled = false;

            for (final String sample : samples) {
                final Iterator<SAMRecordState> iterator = readStates.iterator(sample);
                final List<PileupElement> pile = new ArrayList<PileupElement>(readStates.size(sample));

                int size = 0;                                                           // number of elements in this sample's pileup
                int nDeletions = 0;                                                     // number of deletions in this sample's pileup
                int nMQ0Reads = 0;                                                      // number of MQ0 reads in this sample's pileup (warning: current implementation includes N bases that are MQ0)

                while (iterator.hasNext()) {
                    final SAMRecordState state = iterator.next();                   // state object with the read/offset information
                    final GATKSAMRecord read = (GATKSAMRecord) state.getRead();     // the actual read
                    final CigarOperator op = state.getCurrentCigarOperator();       // current cigar operator
                    final CigarElement nextElement = state.peekForwardOnGenome();   // next cigar element
                    final CigarElement lastElement = state.peekBackwardOnGenome();  // last cigar element
                    final boolean isSingleElementCigar = nextElement == lastElement;
                    final CigarOperator nextOp = nextElement.getOperator();         // next cigar operator
                    final CigarOperator lastOp = lastElement.getOperator();         // last cigar operator
                    int readOffset = state.getReadOffset();                         // the base offset on this read

                    final boolean isBeforeDeletion  = nextOp == CigarOperator.DELETION;
                    final boolean isAfterDeletion   = lastOp == CigarOperator.DELETION;
                    final boolean isBeforeInsertion = nextOp == CigarOperator.INSERTION;
                    final boolean isAfterInsertion  = lastOp == CigarOperator.INSERTION && !isSingleElementCigar;
                    final boolean isNextToSoftClip  = nextOp == CigarOperator.S || (state.getGenomeOffset() == 0 && read.getSoftStart() != read.getAlignmentStart());

                    int nextElementLength = nextElement.getLength();

                    if (op == CigarOperator.N)                                      // N's are never added to any pileup
                        continue;

                    if (op == CigarOperator.D) {
                        // TODO -- LIBS is totally busted for deletions so that reads with Ds right before Is in their CIGAR are broken; must fix
                        if (readInfo.includeReadsWithDeletionAtLoci()) {            // only add deletions to the pileup if we are authorized to do so
                            pile.add(new PileupElement(read, readOffset, true, isBeforeDeletion, isAfterDeletion, isBeforeInsertion, isAfterInsertion, isNextToSoftClip, null, nextOp == CigarOperator.D ? nextElementLength : -1));
                            size++;
                            nDeletions++;
                            if (read.getMappingQuality() == 0)
                                nMQ0Reads++;
                        }
                    }
                    else {
                        if (!filterBaseInRead(read, location.getStart())) {
                            String insertedBaseString = null;
                            if (nextOp == CigarOperator.I) {
                                final int insertionOffset = isSingleElementCigar ? 0 : 1;
                                // TODO -- someone please implement a better fix for the single element insertion CIGAR!
                                if (isSingleElementCigar)
                                    readOffset -= (nextElement.getLength() - 1); // LIBS has passed over the insertion bases!
                                insertedBaseString = new String(Arrays.copyOfRange(read.getReadBases(), readOffset + insertionOffset, readOffset + insertionOffset + nextElement.getLength()));
                            }

                            pile.add(new PileupElement(read, readOffset, false, isBeforeDeletion, isAfterDeletion, isBeforeInsertion, isAfterInsertion, isNextToSoftClip, insertedBaseString, nextElementLength));
                            size++;
                            if (read.getMappingQuality() == 0)
                                nMQ0Reads++;
                        }
                    }
                }

                if (pile.size() != 0)                                             // if this pileup added at least one base, add it to the full pileup
                    fullPileup.put(sample, new ReadBackedPileupImpl(location, pile, size, nDeletions, nMQ0Reads));
            }

            updateReadStates();                                                   // critical - must be called after we get the current state offsets and location
            if (!fullPileup.isEmpty())                                            // if we got reads with non-D/N over the current position, we are done
                nextAlignmentContext = new AlignmentContext(location, new ReadBackedPileupImpl(location, fullPileup), hasBeenSampled);
        }
    }

    // fast testing of position
    private boolean readIsPastCurrentPosition(SAMRecord read) {
        if (readStates.isEmpty())
            return false;
        else {
            SAMRecordState state = readStates.getFirst();
            SAMRecord ourRead = state.getRead();
            return read.getReferenceIndex() > ourRead.getReferenceIndex() || read.getAlignmentStart() > state.getGenomePosition();
        }
    }

    /**
     * Generic place to put per-base filters appropriate to LocusIteratorByState
     *
     * @param rec
     * @param pos
     * @return
     */
    private static boolean filterBaseInRead(GATKSAMRecord rec, long pos) {
        return ReadUtils.isBaseInsideAdaptor(rec, pos);
    }

    private void updateReadStates() {
        for (final String sample : samples) {
            Iterator<SAMRecordState> it = readStates.iterator(sample);
            while (it.hasNext()) {
                SAMRecordState state = it.next();
                CigarOperator op = state.stepForwardOnGenome();
                if (op == null) {
                    // we discard the read only when we are past its end AND indel at the end of the read (if any) was
                    // already processed. Keeping the read state that returned null upon stepForwardOnGenome() is safe
                    // as the next call to stepForwardOnGenome() will return null again AND will clear hadIndel() flag.
                    it.remove();                                                // we've stepped off the end of the object
                }
            }
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    protected class ReadStateManager {
        private final PeekableIterator<SAMRecord> iterator;
        private final SamplePartitioner samplePartitioner;
        private final Map<String, PerSampleReadStateManager> readStatesBySample = new HashMap<String, PerSampleReadStateManager>();
        private int totalReadStates = 0;

        public ReadStateManager(Iterator<SAMRecord> source) {
            this.iterator = new PeekableIterator<SAMRecord>(source);

            for (final String sample : samples) {
                readStatesBySample.put(sample, new PerSampleReadStateManager());
            }

            samplePartitioner = new SamplePartitioner(performDownsampling);
        }

        /**
         * Returns a iterator over all the reads associated with the given sample.  Note that remove() is implemented
         * for this iterator; if present, total read states will be decremented.
         *
         * @param sample The sample.
         * @return Iterator over the reads associated with that sample.
         */
        public Iterator<SAMRecordState> iterator(final String sample) {
            return new Iterator<SAMRecordState>() {
                private Iterator<SAMRecordState> wrappedIterator = readStatesBySample.get(sample).iterator();

                public boolean hasNext() {
                    return wrappedIterator.hasNext();
                }

                public SAMRecordState next() {
                    return wrappedIterator.next();
                }

                public void remove() {
                    wrappedIterator.remove();
                }
            };
        }

        public boolean isEmpty() {
            return totalReadStates == 0;
        }

        /**
         * Retrieves the total number of reads in the manager across all samples.
         *
         * @return Total number of reads over all samples.
         */
        public int size() {
            return totalReadStates;
        }

        /**
         * Retrieves the total number of reads in the manager in the given sample.
         *
         * @param sample The sample.
         * @return Total number of reads in the given sample.
         */
        public int size(final String sample) {
            return readStatesBySample.get(sample).size();
        }

        public SAMRecordState getFirst() {
            for (final String sample : samples) {
                PerSampleReadStateManager reads = readStatesBySample.get(sample);
                if (!reads.isEmpty())
                    return reads.peek();
            }
            return null;
        }

        public boolean hasNext() {
            return totalReadStates > 0 || iterator.hasNext();
        }

        public void collectPendingReads() {
            if (!iterator.hasNext())
                return;

            if (readStates.size() == 0) {
                int firstContigIndex = iterator.peek().getReferenceIndex();
                int firstAlignmentStart = iterator.peek().getAlignmentStart();
                while (iterator.hasNext() && iterator.peek().getReferenceIndex() == firstContigIndex && iterator.peek().getAlignmentStart() == firstAlignmentStart) {
                    samplePartitioner.submitRead(iterator.next());
                }
            } else {
                // Fast fail in the case that the read is past the current position.
                if (readIsPastCurrentPosition(iterator.peek()))
                    return;

                while (iterator.hasNext() && !readIsPastCurrentPosition(iterator.peek())) {
                    samplePartitioner.submitRead(iterator.next());
                }
            }

            samplePartitioner.doneSubmittingReads();

            for (final String sample : samples) {
                Collection<SAMRecord> newReads = samplePartitioner.getReadsForSample(sample);
                PerSampleReadStateManager statesBySample = readStatesBySample.get(sample);
                addReadsToSample(statesBySample, newReads);
            }

            samplePartitioner.reset();
        }

        /**
         * Add reads with the given sample name to the given hanger entry.
         *
         * @param readStates The list of read states to add this collection of reads.
         * @param reads      Reads to add.  Selected reads will be pulled from this source.
         */
        private void addReadsToSample(final PerSampleReadStateManager readStates, final Collection<SAMRecord> reads) {
            if (reads.isEmpty())
                return;

            Collection<SAMRecordState> newReadStates = new LinkedList<SAMRecordState>();

            for (SAMRecord read : reads) {
                SAMRecordState state = new SAMRecordState(read);
                state.stepForwardOnGenome();
                newReadStates.add(state);
            }

            readStates.addStatesAtNextAlignmentStart(newReadStates);
        }

        protected class PerSampleReadStateManager implements Iterable<SAMRecordState> {
            private List<LinkedList<SAMRecordState>> readStatesByAlignmentStart = new LinkedList<LinkedList<SAMRecordState>>();
            private int thisSampleReadStates = 0;
            private Downsampler<LinkedList<SAMRecordState>> levelingDownsampler =
                      performDownsampling ?
                      new LevelingDownsampler<LinkedList<SAMRecordState>, SAMRecordState>(readInfo.getDownsamplingMethod().toCoverage) :
                      null;

            public void addStatesAtNextAlignmentStart(Collection<SAMRecordState> states) {
                if ( states.isEmpty() ) {
                    return;
                }

                readStatesByAlignmentStart.add(new LinkedList<SAMRecordState>(states));
                thisSampleReadStates += states.size();
                totalReadStates += states.size();

                if ( levelingDownsampler != null ) {
                    levelingDownsampler.submit(readStatesByAlignmentStart);
                    levelingDownsampler.signalEndOfInput();

                    thisSampleReadStates -= levelingDownsampler.getNumberOfDiscardedItems();
                    totalReadStates -= levelingDownsampler.getNumberOfDiscardedItems();

                    // use returned List directly rather than make a copy, for efficiency's sake
                    readStatesByAlignmentStart = levelingDownsampler.consumeFinalizedItems();
                    levelingDownsampler.reset();
                }
            }

            public boolean isEmpty() {
                return readStatesByAlignmentStart.isEmpty();
            }

            public SAMRecordState peek() {
                return isEmpty() ? null : readStatesByAlignmentStart.get(0).peek();
            }

            public int size() {
                return thisSampleReadStates;
            }

            public Iterator<SAMRecordState> iterator() {
                return new Iterator<SAMRecordState>() {
                    private Iterator<LinkedList<SAMRecordState>> alignmentStartIterator = readStatesByAlignmentStart.iterator();
                    private LinkedList<SAMRecordState> currentPositionReadStates = null;
                    private Iterator<SAMRecordState> currentPositionReadStatesIterator = null;

                    public boolean hasNext() {
                        return  alignmentStartIterator.hasNext() ||
                                (currentPositionReadStatesIterator != null && currentPositionReadStatesIterator.hasNext());
                    }

                    public SAMRecordState next() {
                        if ( currentPositionReadStatesIterator == null || ! currentPositionReadStatesIterator.hasNext() ) {
                            currentPositionReadStates = alignmentStartIterator.next();
                            currentPositionReadStatesIterator = currentPositionReadStates.iterator();
                        }

                        return currentPositionReadStatesIterator.next();
                    }

                    public void remove() {
                        currentPositionReadStatesIterator.remove();
                        thisSampleReadStates--;
                        totalReadStates--;

                        if ( currentPositionReadStates.isEmpty() ) {
                            alignmentStartIterator.remove();
                        }
                    }
                };
            }
        }
    }

    /**
     * Divides reads by sample and (if requested) does a preliminary downsampling pass with a ReservoirDownsampler.
     *
     * Note: stores reads by sample ID string, not by sample object
     */
    private class SamplePartitioner {
        private Map<String, Downsampler<SAMRecord>> readsBySample;

        public SamplePartitioner( boolean downsampleReads ) {
            readsBySample = new HashMap<String, Downsampler<SAMRecord>>();

            for ( String sample : samples ) {
                readsBySample.put(sample,
                                  downsampleReads ? new ReservoirDownsampler<SAMRecord>(readInfo.getDownsamplingMethod().toCoverage) :
                                                    new PassThroughDownsampler<SAMRecord>());
            }
        }

        public void submitRead(SAMRecord read) {
            String sampleName = read.getReadGroup() != null ? read.getReadGroup().getSample() : null;
            if (readsBySample.containsKey(sampleName))
                readsBySample.get(sampleName).submit(read);
        }

        public void doneSubmittingReads() {
            for ( Map.Entry<String, Downsampler<SAMRecord>> perSampleReads : readsBySample.entrySet() ) {
                perSampleReads.getValue().signalEndOfInput();
            }
        }

        public Collection<SAMRecord> getReadsForSample(String sampleName) {
            if ( ! readsBySample.containsKey(sampleName) )
                throw new NoSuchElementException("Sample name not found");

            return readsBySample.get(sampleName).consumeFinalizedItems();
        }

        public void reset() {
            for ( Map.Entry<String, Downsampler<SAMRecord>> perSampleReads : readsBySample.entrySet() ) {
                perSampleReads.getValue().clear();
                perSampleReads.getValue().reset();
            }
        }
    }
}