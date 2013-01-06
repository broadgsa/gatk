/*
 * Copyright (c) 2012 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.locusiterator;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.downsampling.Downsampler;
import org.broadinstitute.sting.gatk.downsampling.LevelingDownsampler;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * Manages and updates mapping from sample -> List of SAMRecordAlignmentState
 *
 * Optionally can keep track of all of the reads pulled off the iterator and
 * that appeared at any point in the list of SAMRecordAlignmentState for any reads.
 * This functionaly is only possible at this stage, as this object does the popping of
 * reads off the underlying source iterator, and presents only a pileup-like interface
 * of samples -> SAMRecordAlignmentStates.  Reconstructing the unique set of reads
 * used across all pileups is extremely expensive from that data structure.
 *
 * User: depristo
 * Date: 1/5/13
 * Time: 2:02 PM
 */
class ReadStateManager {
    private final List<String> samples;
    private final PeekableIterator<SAMRecord> iterator;
    private final SamplePartitioner samplePartitioner;
    private final Map<String, PerSampleReadStateManager> readStatesBySample = new HashMap<String, PerSampleReadStateManager>();

    private LinkedList<SAMRecord> submittedReads;
    private final boolean keepSubmittedReads;

    private int totalReadStates = 0;

    public ReadStateManager(final Iterator<SAMRecord> source,
                            final List<String> samples,
                            final LIBSDownsamplingInfo LIBSDownsamplingInfo,
                            final boolean keepSubmittedReads) {
        this.samples = samples;
        this.iterator = new PeekableIterator<SAMRecord>(source);

        this.keepSubmittedReads = keepSubmittedReads;
        this.submittedReads = new LinkedList<SAMRecord>();

        for (final String sample : samples) {
            readStatesBySample.put(sample, new PerSampleReadStateManager(LIBSDownsamplingInfo));
        }

        samplePartitioner = new SamplePartitioner(LIBSDownsamplingInfo, samples);
    }

    /**
     * Returns a iterator over all the reads associated with the given sample.  Note that remove() is implemented
     * for this iterator; if present, total read states will be decremented.
     *
     * @param sample The sample.
     * @return Iterator over the reads associated with that sample.
     */
    public Iterator<SAMRecordAlignmentState> iterator(final String sample) {
        return new Iterator<SAMRecordAlignmentState>() {
            private Iterator<SAMRecordAlignmentState> wrappedIterator = readStatesBySample.get(sample).iterator();

            public boolean hasNext() {
                return wrappedIterator.hasNext();
            }

            public SAMRecordAlignmentState next() {
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

    public SAMRecordAlignmentState getFirst() {
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

    // fast testing of position
    private boolean readIsPastCurrentPosition(SAMRecord read) {
        if (isEmpty())
            return false;
        else {
            SAMRecordAlignmentState state = getFirst();
            SAMRecord ourRead = state.getRead();
            return read.getReferenceIndex() > ourRead.getReferenceIndex() || read.getAlignmentStart() > state.getGenomePosition();
        }
    }

    public void collectPendingReads() {
        if (!iterator.hasNext())
            return;

        // the next record in the stream, peeked as to not remove it from the stream
        if ( isEmpty() ) {
            final int firstContigIndex = iterator.peek().getReferenceIndex();
            final int firstAlignmentStart = iterator.peek().getAlignmentStart();
            while (iterator.hasNext() && iterator.peek().getReferenceIndex() == firstContigIndex && iterator.peek().getAlignmentStart() == firstAlignmentStart) {
                submitRead(iterator.next());
            }
        } else {
            // Fast fail in the case that the read is past the current position.
            if (readIsPastCurrentPosition(iterator.peek()))
                return;

            while (iterator.hasNext() && !readIsPastCurrentPosition(iterator.peek())) {
                submitRead(iterator.next());
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
     * Add a read to the sample partitioner, potentially adding it to all submitted reads, if appropriate
     * @param read a non-null read
     */
    @Requires("read != null")
    protected void submitRead(final SAMRecord read) {
        if ( keepSubmittedReads )
            submittedReads.add(read);
        samplePartitioner.submitRead(read);
    }

    /**
     * Transfer current list of submitted reads, clearing old list
     *
     * Takes the maintained list of submitted reads, and transfers it to the caller of this
     * function.  The old list of set to a new, cleanly allocated list so the caller officially
     * owns the list returned by this call.  This is the only way to clear the tracking
     * of submitted reads, if enabled.
     *
     * How to use this function:
     *
     * while ( doing some work unit, such as creating pileup at some locus ):
     *   interact with ReadStateManager in some way to make work unit
     *   readsUsedInPileup = transferSubmittedReads)
     *
     * @throws UnsupportedOperationException if called when keepingSubmittedReads is false
     *
     * @return the current list of submitted reads
     */
    @Ensures({
            "result != null",
            "result != submittedReads" // result and previous submitted reads are not == objects
    })
    public List<SAMRecord> transferSubmittedReads() {
        if ( ! keepSubmittedReads ) throw new UnsupportedOperationException("cannot transferSubmittedReads if you aren't keeping them");

        final List<SAMRecord> prevSubmittedReads = submittedReads;
        this.submittedReads = new LinkedList<SAMRecord>();

        return prevSubmittedReads;
    }

    /**
     * Obtain a pointer to the list of submitted reads.
     *
     * This is not a copy of the list; it is shared with this ReadStateManager.  It should
     * not be modified.  Updates to this ReadStateManager may change the contains of the
     * list entirely.
     *
     * For testing purposes only.
     *
     * Will always be empty if we are are not keepingSubmittedReads
     *
     * @return a non-null list of reads that have been submitted to this ReadStateManager
     */
    @Ensures({"result != null","keepingSubmittedReads || result.isEmpty()"})
    protected List<SAMRecord> getSubmittedReads() {
        return submittedReads;
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

        Collection<SAMRecordAlignmentState> newReadStates = new LinkedList<SAMRecordAlignmentState>();

        for (SAMRecord read : reads) {
            SAMRecordAlignmentState state = new SAMRecordAlignmentState(read);
            state.stepForwardOnGenome();
            newReadStates.add(state);
        }

        readStates.addStatesAtNextAlignmentStart(newReadStates);
    }

    protected class PerSampleReadStateManager implements Iterable<SAMRecordAlignmentState> {
        private List<LinkedList<SAMRecordAlignmentState>> readStatesByAlignmentStart = new LinkedList<LinkedList<SAMRecordAlignmentState>>();
        private final Downsampler<LinkedList<SAMRecordAlignmentState>> levelingDownsampler;

        private int thisSampleReadStates = 0;

        public PerSampleReadStateManager(final LIBSDownsamplingInfo LIBSDownsamplingInfo) {
            this.levelingDownsampler = LIBSDownsamplingInfo.isPerformDownsampling()
                    ? new LevelingDownsampler<LinkedList<SAMRecordAlignmentState>, SAMRecordAlignmentState>(LIBSDownsamplingInfo.getToCoverage())
                    : null;
        }

        public void addStatesAtNextAlignmentStart(Collection<SAMRecordAlignmentState> states) {
            if ( states.isEmpty() ) {
                return;
            }

            readStatesByAlignmentStart.add(new LinkedList<SAMRecordAlignmentState>(states));
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

        public SAMRecordAlignmentState peek() {
            return isEmpty() ? null : readStatesByAlignmentStart.get(0).peek();
        }

        public int size() {
            return thisSampleReadStates;
        }

        public Iterator<SAMRecordAlignmentState> iterator() {
            return new Iterator<SAMRecordAlignmentState>() {
                private Iterator<LinkedList<SAMRecordAlignmentState>> alignmentStartIterator = readStatesByAlignmentStart.iterator();
                private LinkedList<SAMRecordAlignmentState> currentPositionReadStates = null;
                private Iterator<SAMRecordAlignmentState> currentPositionReadStatesIterator = null;

                public boolean hasNext() {
                    return  alignmentStartIterator.hasNext() ||
                            (currentPositionReadStatesIterator != null && currentPositionReadStatesIterator.hasNext());
                }

                public SAMRecordAlignmentState next() {
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
