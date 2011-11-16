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
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.DownsamplingMethod;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.ReservoirDownsampler;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.ExtendedEventPileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileupImpl;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.*;

/** Iterator that traverses a SAM File, accumulating information on a per-locus basis */
public class LocusIteratorByState extends LocusIterator {
    /** our log, which we want to capture anything from this class */
    private static Logger logger = Logger.getLogger(LocusIteratorByState.class);

    // -----------------------------------------------------------------------------------------------------------------
    //
    // member fields
    //
    // -----------------------------------------------------------------------------------------------------------------
    private boolean hasExtendedEvents = false; // will be set to true if at least one read had an indel right before the current position

    /**
     * Used to create new GenomeLocs.
     */
    private final GenomeLocParser genomeLocParser;
    private final ArrayList<String> samples;
    private final ReadStateManager readStates;

    static private class SAMRecordState {
        SAMRecord read;
        int readOffset = -1;     // how far are we offset from the start of the read bases?
        int genomeOffset = -1;   // how far are we offset from the alignment start on the genome?

        Cigar cigar = null;
        int cigarOffset = -1;
        CigarElement curElement = null;
        int nCigarElements = 0;

        // how far are we into a single cigarElement
        int cigarElementCounter = -1;

        // The logical model for generating extended events is as follows: the "record state" implements the traversal
        // along the reference; thus stepForwardOnGenome() returns on every and only on actual reference bases. This
        // can be a (mis)match or a deletion (in the latter case, we still return on every individual reference base the
        // deletion spans). In the extended events mode, the record state also remembers if there was an insertion, or
        // if the deletion just started *right before* the current reference base the record state is pointing to upon the return from
        // stepForwardOnGenome(). The next call to stepForwardOnGenome() will clear that memory (as we remember only extended
        // events immediately preceding the current reference base).

        boolean generateExtendedEvents = true; // should we generate an additional,  special pile for indels between the ref bases?
        // the only purpose of this flag is to shield away a few additional lines of code
        // when extended piles are not needed, it may not be even worth it...
        byte[] insertedBases = null; // remember full inserted sequence if we are generating piles of extended events (indels)
        int eventLength = -1; // will be set to the length of insertion/deletion if we are generating piles of extended events
        byte eventDelayedFlag = 0; // will be set to non-0 if there was an event (indel) right before the
        // current base on the ref. We use a counter-like variable here since clearing the indel event is
        // delayed by one base, so we need to remember how long ago we have seen the actual event
        int eventStart = -1; // where on the read the extended event starts (i.e. the last position on the read prior to the
        // event, or -1 if alignment starts with an insertion); this one is easy to recompute on the fly,
        // we cache it here mainly for convenience


        public SAMRecordState(SAMRecord read, boolean extended) {
            this.read = read;
            cigar = read.getCigar();
            nCigarElements = cigar.numCigarElements();
            generateExtendedEvents = extended;

            //System.out.printf("Creating a SAMRecordState: %s%n", this);
        }

        public SAMRecord getRead() { return read; }

        /**
         * What is our current offset in the read's bases that aligns us with the reference genome?
         *
         * @return
         */
        public int getReadOffset() { return readOffset; }

        /**
         * What is the current offset w.r.t. the alignment state that aligns us to the readOffset?
         *
         * @return
         */
        public int getGenomeOffset() { return genomeOffset; }

        public int getGenomePosition() { return read.getAlignmentStart() + getGenomeOffset(); }

        public GenomeLoc getLocation(GenomeLocParser genomeLocParser) {
            return genomeLocParser.createGenomeLoc(read.getReferenceName(), getGenomePosition());
        }

        public CigarOperator getCurrentCigarOperator() {
            return curElement.getOperator();
        }

        /** Returns true if we just stepped over insertion/into a deletion prior to the last return from stepForwardOnGenome.
         *
         * @return
         */
        public boolean hadIndel() {
            return ( eventLength > 0 );
        }

        public int getEventLength() { return eventLength; }

        public byte[] getEventBases() { return insertedBases; }

        public int getReadEventStartOffset() { return eventStart; }

        public String toString() {
            return String.format("%s ro=%d go=%d co=%d cec=%d %s", read.getReadName(), readOffset, genomeOffset, cigarOffset, cigarElementCounter, curElement);
        }

        public CigarOperator stepForwardOnGenome() {
            // we enter this method with readOffset = index of the last processed base on the read
            // (-1 if we did not process a single base yet); this can be last matching base, or last base of an insertion


            if ( curElement == null || ++cigarElementCounter > curElement.getLength() ) {
                cigarOffset++;
                if ( cigarOffset < nCigarElements ) {
                    curElement = cigar.getCigarElement(cigarOffset);
                    cigarElementCounter = 0;
                    // next line: guards against cigar elements of length 0; when new cigar element is retrieved,
                    // we reenter in order to re-check cigarElementCounter against curElement's length
                    return stepForwardOnGenome();
                } else {
                    // Reads that contain indels model the genomeOffset as the following base in the reference.  Because
                    // we fall into this else block only when indels end the read, increment genomeOffset  such that the
                    // current offset of this read is the next ref base after the end of the indel.  This position will
                    // model a point on the reference somewhere after the end of the read.
                    genomeOffset++; // extended events need that. Logically, it's legal to advance the genomic offset here:
                                    // we do step forward on the ref, and by returning null we also indicate that we are past the read end.

                    if ( generateExtendedEvents && eventDelayedFlag > 0 ) {

                        // if we had an indel right before the read ended (i.e. insertion was the last cigar element),
                        // we keep it until next reference base; then we discard it and this will allow the LocusIterator to
                        // finally discard this read
                        eventDelayedFlag--;
                        if ( eventDelayedFlag == 0 )  {
                            eventLength = -1; // reset event when we are past it
                            insertedBases = null;
                            eventStart = -1;
                        }
                    }

                    return null;
                }
            }


            boolean done = false;
            switch (curElement.getOperator()) {
                case H : // ignore hard clips
                case P : // ignore pads
                    cigarElementCounter = curElement.getLength();
                    break;
                case I : // insertion w.r.t. the reference
                    if ( generateExtendedEvents ) {
                        // we see insertions only once, when we step right onto them; the position on the read is scrolled
                        // past the insertion right after that
                        if ( eventDelayedFlag > 1 ) throw new UserException.MalformedBAM(read, "Adjacent I/D events in read "+read.getReadName());
                        insertedBases = Arrays.copyOfRange(read.getReadBases(),readOffset+1,readOffset+1+curElement.getLength());
                        eventLength = curElement.getLength() ;
                        eventStart = readOffset;
                        eventDelayedFlag = 2; // insertion causes re-entry into stepForwardOnGenome, so we set the delay to 2
//                        System.out.println("Inserted "+(new String (insertedBases)) +" after "+readOffset);
                    } // continue onto the 'S' case !
                case S : // soft clip
                    cigarElementCounter = curElement.getLength();
                    readOffset += curElement.getLength();
                    break;
                case D : // deletion w.r.t. the reference
                    if ( generateExtendedEvents ) {
                        if ( cigarElementCounter == 1) {
                            // generate an extended event only if we just stepped into the deletion (i.e. don't
                            // generate the event at every deleted position on the ref, that's what cigarElementCounter==1 is for!)
                            if ( eventDelayedFlag > 1 ) throw new UserException.MalformedBAM(read, "Adjacent I/D events in read "+read.getReadName());
                            eventLength = curElement.getLength();
                            eventDelayedFlag = 2; // deletion on the ref causes an immediate return, so we have to delay by 1 only
                            eventStart = readOffset;
                            insertedBases = null;
//                            System.out.println("Deleted "+eventLength +" bases after "+readOffset);
                        }
                    }
                    // should be the same as N case
                    genomeOffset++;
                    done = true;
                    break;
                case N : // reference skip (looks and gets processed just like a "deletion", just different logical meaning)
                    genomeOffset++;
                    done = true;
                    break;
                case M :
                    readOffset++;
                    genomeOffset++;
                    done = true;
                    break;
                default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + curElement.getOperator());
            }

            if ( generateExtendedEvents ) {
                if ( eventDelayedFlag > 0 && done ) {
                    // if we did make a successful step on the ref, decrement delayed flag. If, upon the decrementthe,
                    // the flag is 1, we are standing on the reference base right after the indel (so we have to keep it).
                    // Otherwise, we are away from the previous indel and have to clear our memories...
                    eventDelayedFlag--; // when we notice an indel, we set delayed flag to 2, so now
                                    // if eventDelayedFlag == 1, an indel occured right before the current base
                    if ( eventDelayedFlag == 0 ) {
                        eventLength = -1; // reset event when we are past it
                        insertedBases = null;
                        eventStart = -1;
                    }
                }
            }

            return done ? curElement.getOperator() : stepForwardOnGenome();
        }
    }

    //final boolean DEBUG = false;
    //final boolean DEBUG2 = false && DEBUG;
    private ReadProperties readInfo;
    private AlignmentContext nextAlignmentContext;

    // -----------------------------------------------------------------------------------------------------------------
    //
    // constructors and other basic operations
    //
    // -----------------------------------------------------------------------------------------------------------------

    public LocusIteratorByState(final Iterator<SAMRecord> samIterator, ReadProperties readInformation, GenomeLocParser genomeLocParser, Collection<String> samples ) {
        this.readInfo = readInformation;
        this.genomeLocParser = genomeLocParser;
        this.samples = new ArrayList<String>(samples);
        this.readStates = new ReadStateManager(samIterator,readInformation.getDownsamplingMethod());

        // currently the GATK expects this LocusIteratorByState to accept empty sample lists, when
        // there's no read data.  So we need to throw this error only when samIterator.hasNext() is true
        if ( this.samples.isEmpty() && samIterator.hasNext() ) {
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
        if(!hasNext())
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
        while(nextAlignmentContext == null && readStates.hasNext()) {
            // this call will set hasExtendedEvents to true if it picks up a read with indel right before the current position on the ref:
            readStates.collectPendingReads();

            int size = 0;
            int nDeletions = 0;
            int nInsertions = 0;
            int nMQ0Reads = 0;


            // if extended events are requested, and if previous traversal step brought us over an indel in
            // at least one read, we emit extended pileup (making sure that it is associated with the previous base,
            // i.e. the one right *before* the indel) and do NOT shift the current position on the ref.
            // In this case, the subsequent call to next() will emit the normal pileup at the current base
            // and shift the position.
            if (readInfo.generateExtendedEvents() && hasExtendedEvents) {
                Map<String,ReadBackedExtendedEventPileupImpl> fullExtendedEventPileup = new HashMap<String,ReadBackedExtendedEventPileupImpl>();

                // get current location on the reference and decrement it by 1: the indels we just stepped over
                // are associated with the *previous* reference base
                GenomeLoc loc = genomeLocParser.incPos(getLocation(),-1);

                boolean hasBeenSampled = false;
                for(final String sample: samples) {
                    Iterator<SAMRecordState> iterator = readStates.iterator(sample);
                    List<ExtendedEventPileupElement> indelPile = new ArrayList<ExtendedEventPileupElement>(readStates.size(sample));
                    hasBeenSampled |= loc.getStart() <= readStates.getDownsamplingExtent(sample);

                    size = 0;
                    nDeletions = 0;
                    nInsertions = 0;
                    nMQ0Reads = 0;
                    int maxDeletionLength = 0;

                    while(iterator.hasNext()) {
                        SAMRecordState state = iterator.next();
                        if ( state.hadIndel() ) {
                            size++;
                            if ( state.getEventBases() == null ) {
                                nDeletions++;
                                maxDeletionLength = Math.max(maxDeletionLength,state.getEventLength());
                            }
                            else nInsertions++;
                            indelPile.add ( new ExtendedEventPileupElement((GATKSAMRecord) state.getRead(), state.getReadEventStartOffset(), state.getEventLength(), state.getEventBases()) );

                        }   else {
                            // HACK: The readahead mechanism for LocusIteratorByState will effectively read past the current position
                            //       and add in extra reads that start after this indel.  Skip these reads.
                            //       My belief at this moment after empirically looking at read->ref alignment is that, in a cigar string
                            //       like 1I76M, the first insertion is between alignment start-1 and alignment start, so we shouldn't be
                            //       filtering these out.
                            // TODO: UPDATE!  Eric tells me that we *might* want reads adjacent to the pileup in the pileup.  Strike this block.
                            //if(state.getRead().getAlignmentStart() > loc.getStart())
                            //    continue;

                            if ( state.getCurrentCigarOperator() != CigarOperator.N ) {
                                // this read has no indel associated with the previous position on the ref;
                                // we count this read in only if it has actual bases, not N span...
                                if ( state.getCurrentCigarOperator() != CigarOperator.D || readInfo.includeReadsWithDeletionAtLoci() ) {

                                    // if cigar operator is D but the read has no extended event reported (that's why we ended
                                    // up in this branch), it means that we are currently inside a deletion that started earlier;
                                    // we count such reads (with a longer deletion spanning over a deletion at the previous base we are
                                    // about to report) only if includeReadsWithDeletionAtLoci is true.
                                    size++;
                                    indelPile.add ( new ExtendedEventPileupElement((GATKSAMRecord) state.getRead(), state.getReadOffset()-1, -1) // length=-1 --> noevent
                                            );
                                }
                            }
                        }
                        if ( state.getRead().getMappingQuality() == 0 ) {
                            nMQ0Reads++;
                        }
                    }
                    if( indelPile.size() != 0 ) fullExtendedEventPileup.put(sample,new ReadBackedExtendedEventPileupImpl(loc,indelPile,size,maxDeletionLength,nInsertions,nDeletions,nMQ0Reads));
                }
                hasExtendedEvents = false; // we are done with extended events prior to current ref base
//                System.out.println("Indel(s) at "+loc);
//               for ( ExtendedEventPileupElement pe : indelPile ) { if ( pe.isIndel() ) System.out.println("  "+pe.toString()); }
                nextAlignmentContext = new AlignmentContext(loc, new ReadBackedExtendedEventPileupImpl(loc, fullExtendedEventPileup), hasBeenSampled);
            }  else {
                GenomeLoc location = getLocation();
                Map<String,ReadBackedPileupImpl> fullPileup = new HashMap<String,ReadBackedPileupImpl>();

                boolean hasBeenSampled = false;
                for(final String sample: samples) {
                    Iterator<SAMRecordState> iterator = readStates.iterator(sample);
                    List<PileupElement> pile = new ArrayList<PileupElement>(readStates.size(sample));
                    hasBeenSampled |= location.getStart() <= readStates.getDownsamplingExtent(sample);

                    size = 0;
                    nDeletions = 0;
                    nMQ0Reads = 0;

                    while(iterator.hasNext()) {
                        SAMRecordState state = iterator.next();
                        if ( state.getCurrentCigarOperator() != CigarOperator.D && state.getCurrentCigarOperator() != CigarOperator.N ) {
                            if ( filterBaseInRead(state.getRead(), location.getStart()) ) {
                                //discarded_bases++;
                                //printStatus("Adaptor bases", discarded_adaptor_bases);
                                continue;
                            } else {
                                //observed_bases++;
                                pile.add(new PileupElement((GATKSAMRecord) state.getRead(), state.getReadOffset()));
                                size++;
                            }
                        } else if ( readInfo.includeReadsWithDeletionAtLoci() && state.getCurrentCigarOperator() != CigarOperator.N ) {
                            size++;
                            pile.add(new PileupElement((GATKSAMRecord) state.getRead(), -1));
                            nDeletions++;
                        }

                        if ( state.getRead().getMappingQuality() == 0 ) {
                            nMQ0Reads++;
                        }
                    }

                    if( pile.size() != 0 )
                        fullPileup.put(sample,new ReadBackedPileupImpl(location,pile,size,nDeletions,nMQ0Reads));
                }

                updateReadStates(); // critical - must be called after we get the current state offsets and location
                // if we got reads with non-D/N over the current position, we are done
                if ( !fullPileup.isEmpty() ) nextAlignmentContext = new AlignmentContext(location, new ReadBackedPileupImpl(location,fullPileup),hasBeenSampled);
            }
        }
    }

    // fast testing of position
    private boolean readIsPastCurrentPosition(SAMRecord read) {
        if ( readStates.isEmpty() )
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
    private static boolean filterBaseInRead(SAMRecord rec, long pos) {
        return ReadUtils.readPairBaseOverlapType(rec, pos) == ReadUtils.OverlapType.IN_ADAPTOR;
    }

    private void updateReadStates() {
        for(final String sample: samples) {
            Iterator<SAMRecordState> it = readStates.iterator(sample);
            while ( it.hasNext() ) {
                SAMRecordState state = it.next();
                CigarOperator op = state.stepForwardOnGenome();
                if ( state.hadIndel() && readInfo.generateExtendedEvents() ) hasExtendedEvents = true;
                else {
                    // we discard the read only when we are past its end AND indel at the end of the read (if any) was
                    // already processed. Keeping the read state that retunred null upon stepForwardOnGenome() is safe
                    // as the next call to stepForwardOnGenome() will return null again AND will clear hadIndel() flag.
                    if ( op == null ) { // we've stepped off the end of the object
                        //if (DEBUG) logger.debug(String.format("   removing read %s at %d", state.getRead().getReadName(), state.getRead().getAlignmentStart()));
                        it.remove();
                    }
                }
            }
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    private class ReadStateManager  {
        private final PeekableIterator<SAMRecord> iterator;
        private final DownsamplingMethod downsamplingMethod;
        private final SamplePartitioner samplePartitioner;
        private final Map<String,PerSampleReadStateManager> readStatesBySample = new HashMap<String,PerSampleReadStateManager>();
        private final int targetCoverage;
        private int totalReadStates = 0;

        public ReadStateManager(Iterator<SAMRecord> source, DownsamplingMethod downsamplingMethod) {
            this.iterator = new PeekableIterator<SAMRecord>(source);
            this.downsamplingMethod = downsamplingMethod.type != null ? downsamplingMethod : DownsamplingMethod.NONE;
            switch(this.downsamplingMethod.type) {
                case BY_SAMPLE:
                    if(downsamplingMethod.toCoverage == null)
                        throw new UserException.BadArgumentValue("dcov", "Downsampling coverage (-dcov) must be specified when downsampling by sample");
                    this.targetCoverage = downsamplingMethod.toCoverage;
                    break;
                default:
                    this.targetCoverage = Integer.MAX_VALUE;
            }

            Map<String,ReadSelector> readSelectors = new HashMap<String,ReadSelector>();
            for(final String sample: samples) {
                readStatesBySample.put(sample,new PerSampleReadStateManager());
                readSelectors.put(sample,downsamplingMethod.type == DownsampleType.BY_SAMPLE ? new NRandomReadSelector(null,targetCoverage) : new AllReadsSelector());
            }

            samplePartitioner = new SamplePartitioner(readSelectors);
        }

        /**
         * Returns a iterator over all the reads associated with the given sample.  Note that remove() is implemented
         * for this iterator; if present, total read states will be decremented.
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
                    totalReadStates--;
                }
            };
        }

        public boolean isEmpty() {
            return totalReadStates == 0;
        }

        /**
         * Retrieves the total number of reads in the manager across all samples.
         * @return Total number of reads over all samples.
         */
        public int size() {
            return totalReadStates;
        }

        /**
         * Retrieves the total number of reads in the manager in the given sample.
         * @param sample The sample.
         * @return Total number of reads in the given sample.
         */
        public int size(final String sample) {
            return readStatesBySample.get(sample).size();
        }

        /**
         * The extent of downsampling; basically, the furthest base out which has 'fallen
         * victim' to the downsampler.
         * @param sample Sample, downsampled independently.
         * @return Integer stop of the furthest undownsampled region.
         */
        public int getDownsamplingExtent(final String sample) {
            return readStatesBySample.get(sample).getDownsamplingExtent();
        }

        public SAMRecordState getFirst() {
            for(final String sample: samples) {
                PerSampleReadStateManager reads = readStatesBySample.get(sample);
                if(!reads.isEmpty())
                    return reads.peek();
            }
            return null;
        }

        public boolean hasNext() {
            return totalReadStates > 0 || iterator.hasNext();
        }

        public void collectPendingReads() {
            if(!iterator.hasNext())
                return;

            if(readStates.size() == 0) {
                int firstContigIndex = iterator.peek().getReferenceIndex();
                int firstAlignmentStart = iterator.peek().getAlignmentStart();
                while(iterator.hasNext() && iterator.peek().getReferenceIndex() == firstContigIndex && iterator.peek().getAlignmentStart() == firstAlignmentStart) {
                    samplePartitioner.submitRead(iterator.next());
                }
            }
            else {
                // Fast fail in the case that the read is past the current position.
                if(readIsPastCurrentPosition(iterator.peek()))
                    return;

                while (iterator.hasNext() && !readIsPastCurrentPosition(iterator.peek())) {
                    samplePartitioner.submitRead(iterator.next());
                }
            }
            samplePartitioner.complete();

            for(final String sample: samples) {
                ReadSelector aggregator = samplePartitioner.getSelectedReads(sample);

                Collection<SAMRecord> newReads = new ArrayList<SAMRecord>(aggregator.getSelectedReads());

                PerSampleReadStateManager statesBySample = readStatesBySample.get(sample);
                int numReads = statesBySample.size();
                int downsamplingExtent = aggregator.getDownsamplingExtent();

                if(numReads+newReads.size()<=targetCoverage || downsamplingMethod.type==DownsampleType.NONE) {
                    long readLimit = aggregator.getNumReadsSeen();
                    addReadsToSample(statesBySample,newReads,readLimit);
                    statesBySample.specifyNewDownsamplingExtent(downsamplingExtent);
                }
                else {
                    int[] counts = statesBySample.getCountsPerAlignmentStart();
                    int[] updatedCounts = new int[counts.length];
                    System.arraycopy(counts,0,updatedCounts,0,counts.length);

                    boolean readPruned = true;
                    while(numReads+newReads.size()>targetCoverage && readPruned) {
                        readPruned = false;
                        for(int alignmentStart=updatedCounts.length-1;numReads+newReads.size()>targetCoverage&&alignmentStart>=0;alignmentStart--) {
                            if(updatedCounts[alignmentStart] > 1) {
                                updatedCounts[alignmentStart]--;
                                numReads--;
                                readPruned = true;
                            }
                        }
                    }

                    if(numReads == targetCoverage) {
                        updatedCounts[0]--;
                        numReads--;
                    }

                    BitSet toPurge = new BitSet(readStates.size());
                    int readOffset = 0;

                    for(int i = 0; i < updatedCounts.length; i++) {
                        int n = counts[i];
                        int k = updatedCounts[i];

                        for(Integer purgedElement: MathUtils.sampleIndicesWithoutReplacement(n,n-k))
                            toPurge.set(readOffset+purgedElement);

                        readOffset += counts[i];
                    }
                    downsamplingExtent = Math.max(downsamplingExtent,statesBySample.purge(toPurge));
                    
                    addReadsToSample(statesBySample,newReads,targetCoverage-numReads);
                    statesBySample.specifyNewDownsamplingExtent(downsamplingExtent);
                }
            }
            samplePartitioner.reset();
        }

        /**
         * Add reads with the given sample name to the given hanger entry.
         * @param readStates The list of read states to add this collection of reads.
         * @param reads Reads to add.  Selected reads will be pulled from this source.
         * @param maxReads Maximum number of reads to add.
         */
        private void addReadsToSample(final PerSampleReadStateManager readStates, final Collection<SAMRecord> reads, final long maxReads) {
            if(reads.isEmpty())
                return;

            Collection<SAMRecordState> newReadStates = new LinkedList<SAMRecordState>();
            int readCount = 0;
            for(SAMRecord read: reads) {
                if(readCount < maxReads) {
                    SAMRecordState state = new SAMRecordState(read, readInfo.generateExtendedEvents());
                    state.stepForwardOnGenome();
                    newReadStates.add(state);
                    // TODO: What if we downsample the extended events away?
                    if (state.hadIndel()) hasExtendedEvents = true;
                    readCount++;
                }
            }
            readStates.addStatesAtNextAlignmentStart(newReadStates);
        }

        private class PerSampleReadStateManager implements Iterable<SAMRecordState> {
            private final Queue<SAMRecordState> readStates = new LinkedList<SAMRecordState>();
            private final Deque<Counter> readStateCounter = new LinkedList<Counter>();
            private int downsamplingExtent = 0;

            public void addStatesAtNextAlignmentStart(Collection<SAMRecordState> states) {
                readStates.addAll(states);
                readStateCounter.add(new Counter(states.size()));
                totalReadStates += states.size();
            }

            public boolean isEmpty() {
                return readStates.isEmpty();
            }

            public SAMRecordState peek() {
                return readStates.peek();
            }

            public int size() {
                return readStates.size();
            }

            public void specifyNewDownsamplingExtent(int downsamplingExtent) {
                this.downsamplingExtent = Math.max(this.downsamplingExtent,downsamplingExtent);
            }

            public int getDownsamplingExtent() {
                return downsamplingExtent;
            }

            public int[] getCountsPerAlignmentStart() {
                int[] counts = new int[readStateCounter.size()];
                int index = 0;
                for(Counter counter: readStateCounter)
                    counts[index++] = counter.getCount();
                return counts;
            }

            public Iterator<SAMRecordState> iterator() {
                return new Iterator<SAMRecordState>() {
                    private Iterator<SAMRecordState> wrappedIterator = readStates.iterator();

                    public boolean hasNext() {
                        return wrappedIterator.hasNext();
                    }

                    public SAMRecordState next() {
                        return wrappedIterator.next();
                    }

                    public void remove() {
                        wrappedIterator.remove();
                        Counter counter = readStateCounter.peek();
                        counter.decrement();
                        if(counter.getCount() == 0)
                            readStateCounter.remove();
                    }
                };
            }

            /**
             * Purge the given elements from the bitset.  If an element in the bitset is true, purge
             * the corresponding read state.
             * @param elements bits from the set to purge.
             * @return the extent of the final downsampled read.
             */
            public int purge(final BitSet elements) {
                int downsamplingExtent = 0;

                if(elements.isEmpty() || readStates.isEmpty()) return downsamplingExtent;

                Iterator<SAMRecordState> readStateIterator = readStates.iterator();

                Iterator<Counter> counterIterator = readStateCounter.iterator();
                Counter currentCounter = counterIterator.next();

                int readIndex = 0;
                long alignmentStartCounter = currentCounter.getCount();

                int toPurge = elements.nextSetBit(0);
                int removedCount = 0;

                while(readStateIterator.hasNext() && toPurge >= 0) {
                    SAMRecordState state = readStateIterator.next();
                    downsamplingExtent = Math.max(downsamplingExtent,state.getRead().getAlignmentEnd());

                    if(readIndex == toPurge) {
                        readStateIterator.remove();
                        currentCounter.decrement();
                        if(currentCounter.getCount() == 0)
                            counterIterator.remove();
                        removedCount++;
                        toPurge = elements.nextSetBit(toPurge+1);
                    }

                    readIndex++;
                    alignmentStartCounter--;
                    if(alignmentStartCounter == 0 && counterIterator.hasNext()) {
                        currentCounter = counterIterator.next();
                        alignmentStartCounter = currentCounter.getCount();
                    }
                }

                totalReadStates -= removedCount;

                return downsamplingExtent;
            }
        }
    }

    /**
     * Note: assuming that, whenever we downsample, we downsample to an integer capacity.
     */
    static private class Counter {
        private int count;

        public Counter(int count) {
            this.count = count;
        }

        public int getCount() {
            return count;
        }

        public void decrement() {
            count--;
        }
    }
}

/**
 * Selects reads passed to it based on a criteria decided through inheritance.
 * TODO: This is a temporary abstraction until we can get rid of this downsampling implementation and the mrl option.  Get rid of this.
 */
interface ReadSelector {
    /**
     * All previous selectors in the chain have allowed this read.  Submit it to this selector for consideration.
     * @param read the read to evaluate.
     */
    public void submitRead(SAMRecord read);

    /**
     * A previous selector has deemed this read unfit.  Notify this selector so that this selector's counts are valid.
     * @param read the read previously rejected.
     */
    public void notifyReadRejected(SAMRecord read);

    /**
     * Signal the selector that read additions are complete.
     */
    public void complete();

    /**
     * Retrieve the number of reads seen by this selector so far.
     * @return number of reads seen.
     */
    public long getNumReadsSeen();

    /**
     * Return the number of reads accepted by this selector so far.
     * @return number of reads selected.
     */
    public long getNumReadsSelected();

    /**
     * Gets the locus at which the last of the downsampled reads selected by this selector ends.  The value returned will be the
     * last aligned position from this selection to which a downsampled read aligns -- in other words, if a read is thrown out at
     * position 3 whose cigar string is 76M, the value of this parameter will be 78.
     * @return If any read has been downsampled, this will return the last aligned base of the longest alignment.  Else, 0.
     */
    public int getDownsamplingExtent();

    /**
     * Get the reads selected by this selector.
     * @return collection of reads selected by this selector.
     */
    public Collection<SAMRecord> getSelectedReads();

    /**
     * Reset this collection to its pre-gathered state.
     */
    public void reset();
}

/**
 * Select every read passed in.
 */
class AllReadsSelector implements ReadSelector {
    private Collection<SAMRecord> reads = new LinkedList<SAMRecord>();
    private long readsSeen = 0;
    private int downsamplingExtent = 0;

    public void submitRead(SAMRecord read) {
        reads.add(read);
        readsSeen++;
    }

    public void notifyReadRejected(SAMRecord read) {
        readsSeen++;
        downsamplingExtent = Math.max(downsamplingExtent,read.getAlignmentEnd());
    }

    public void complete() {
        // NO-OP.
    }

    public long getNumReadsSeen() {
        return readsSeen;
    }

    public long getNumReadsSelected() {
        return readsSeen;
    }

    public int getDownsamplingExtent() {
        return downsamplingExtent;
    }

    public Collection<SAMRecord> getSelectedReads() {
        return reads;
    }

    public void reset() {
        reads.clear();
        readsSeen = 0;
        downsamplingExtent = 0;
    }
}


/**
 * Select N reads randomly from the input stream.
 */
class NRandomReadSelector implements ReadSelector {
    private final ReservoirDownsampler<SAMRecord> reservoir;
    private final ReadSelector chainedSelector;
    private long readsSeen = 0;
    private int downsamplingExtent = 0;    

    public NRandomReadSelector(ReadSelector chainedSelector, long readLimit) {
        this.reservoir = new ReservoirDownsampler<SAMRecord>((int)readLimit);
        this.chainedSelector = chainedSelector;
    }

    public void submitRead(SAMRecord read) {
        SAMRecord displaced = reservoir.add(read);
        if(displaced != null && chainedSelector != null) {
            chainedSelector.notifyReadRejected(read);
            downsamplingExtent = Math.max(downsamplingExtent,read.getAlignmentEnd());
        }
        readsSeen++;
    }

    public void notifyReadRejected(SAMRecord read) {
        readsSeen++;
    }

    public void complete() {
        for(SAMRecord read: reservoir.getDownsampledContents())
            chainedSelector.submitRead(read);
        if(chainedSelector != null)
            chainedSelector.complete();
    }


    public long getNumReadsSeen() {
        return readsSeen;
    }

    public long getNumReadsSelected() {
        return reservoir.size();
    }

    public int getDownsamplingExtent() {
        return downsamplingExtent;
    }    

    public Collection<SAMRecord> getSelectedReads() {
        return reservoir.getDownsampledContents();
    }

    public void reset() {
        reservoir.clear();
        downsamplingExtent = 0;
        if(chainedSelector != null)
            chainedSelector.reset();
    }
}

/**
 * Note: stores reads by sample ID string, not by sample object
 */
class SamplePartitioner implements ReadSelector {
    private final Map<String,ReadSelector> readsBySample;
    private long readsSeen = 0;

    public SamplePartitioner(Map<String,ReadSelector> readSelectors) {
        readsBySample = readSelectors;
    }

    public void submitRead(SAMRecord read) {
        String sampleName = read.getReadGroup()!=null ? read.getReadGroup().getSample() : null;
        if(readsBySample.containsKey(sampleName))
            readsBySample.get(sampleName).submitRead(read);
        readsSeen++;
    }

    public void notifyReadRejected(SAMRecord read) {
        String sampleName = read.getReadGroup()!=null ? read.getReadGroup().getSample() : null;
        if(readsBySample.containsKey(sampleName))
            readsBySample.get(sampleName).notifyReadRejected(read);
        readsSeen++;
    }

    public void complete() {
        // NO-OP.
    }

    public long getNumReadsSeen() {
        return readsSeen;
    }

    public long getNumReadsSelected() {
        return readsSeen;
    }

    public int getDownsamplingExtent() {
        int downsamplingExtent = 0;
        for(ReadSelector storage: readsBySample.values())
            downsamplingExtent = Math.max(downsamplingExtent,storage.getDownsamplingExtent());
        return downsamplingExtent;
    }
    
    public Collection<SAMRecord> getSelectedReads() {
        throw new UnsupportedOperationException("Cannot directly get selected reads from a read partitioner.");
    }

    public ReadSelector getSelectedReads(String sampleName) {
        if(!readsBySample.containsKey(sampleName))
            throw new NoSuchElementException("Sample name not found");
        return readsBySample.get(sampleName);
    }

    public void reset() {
        for(ReadSelector storage: readsBySample.values())
            storage.reset();
        readsSeen = 0;
    }

}
