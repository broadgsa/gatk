/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.downsampling;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;

import java.util.*;


/**
 * GATKSAMIterator wrapper around our generic reads downsampler interface
 * that downsamples reads for each sample independently, and then re-assembles
 * the reads back into a single merged stream.
 *
 * @author David Roazen
 */
public class PerSampleDownsamplingReadsIterator implements GATKSAMIterator {

    private GATKSAMIterator nestedSAMIterator;
    private ReadsDownsamplerFactory<SAMRecord> downsamplerFactory;
    private Map<String, ReadsDownsampler<SAMRecord>> perSampleDownsamplers;
    private PriorityQueue<SAMRecord> orderedDownsampledReadsCache;
    private SAMRecord nextRead = null;
    private SAMRecordComparator readComparator = new SAMRecordCoordinateComparator();
    private SAMRecord earliestPendingRead = null;
    private ReadsDownsampler<SAMRecord> earliestPendingDownsampler = null;

    // Initial size of our cache of finalized reads
    private static final int DOWNSAMPLED_READS_INITIAL_CACHE_SIZE = 4096;

    // The number of positional changes that can occur in the read stream before all downsamplers
    // should be informed of the current position (guards against samples with relatively sparse reads
    // getting stuck in a pending state):
    private static final int DOWNSAMPLER_POSITIONAL_UPDATE_INTERVAL = 3;   // TODO: experiment with this value

    /**
     * @param iter wrapped iterator from which this iterator will pull reads
     * @param downsamplerFactory factory used to create new downsamplers as needed
     */
    public PerSampleDownsamplingReadsIterator( GATKSAMIterator iter, ReadsDownsamplerFactory<SAMRecord> downsamplerFactory ) {
        nestedSAMIterator = iter;
        this.downsamplerFactory = downsamplerFactory;
        perSampleDownsamplers = new HashMap<String, ReadsDownsampler<SAMRecord>>();
        orderedDownsampledReadsCache = new PriorityQueue<SAMRecord>(DOWNSAMPLED_READS_INITIAL_CACHE_SIZE, readComparator);

        advanceToNextRead();
    }

    public boolean hasNext() {
        return nextRead != null;
    }

    public SAMRecord next() {
        if ( nextRead == null ) {
            throw new NoSuchElementException("next() called when there are no more items");
        }

        SAMRecord toReturn = nextRead;
        advanceToNextRead();

        return toReturn;
    }

    private void advanceToNextRead() {
        if ( ! readyToReleaseReads() && ! fillDownsampledReadsCache() ) {
            nextRead = null;
        }
        else {
            nextRead = orderedDownsampledReadsCache.poll();
        }
    }

    private boolean readyToReleaseReads() {
        if ( orderedDownsampledReadsCache.isEmpty() ) {
            return false;
        }

        return earliestPendingRead == null ||
               readComparator.compare(orderedDownsampledReadsCache.peek(), earliestPendingRead) <= 0;
    }

    private boolean fillDownsampledReadsCache() {
        SAMRecord prevRead = null;
        int numPositionalChanges = 0;

        // Continue submitting reads to the per-sample downsamplers until the read at the top of the priority queue
        // can be released without violating global sort order
        while ( nestedSAMIterator.hasNext() && ! readyToReleaseReads() ) {
            SAMRecord read = nestedSAMIterator.next();
            String sampleName = read.getReadGroup() != null ? read.getReadGroup().getSample() : null;

            ReadsDownsampler<SAMRecord> thisSampleDownsampler = perSampleDownsamplers.get(sampleName);
            if ( thisSampleDownsampler == null ) {
                thisSampleDownsampler = downsamplerFactory.newInstance();
                perSampleDownsamplers.put(sampleName, thisSampleDownsampler);
            }

            thisSampleDownsampler.submit(read);
            processFinalizedAndPendingItems(thisSampleDownsampler);

            if ( prevRead != null && prevRead.getAlignmentStart() != read.getAlignmentStart() ) {
                numPositionalChanges++;
            }

            // Periodically inform all downsamplers of the current position in the read stream. This is
            // to prevent downsamplers for samples with sparser reads than others from getting stuck too
            // long in a pending state.
            if ( numPositionalChanges > 0 && numPositionalChanges % DOWNSAMPLER_POSITIONAL_UPDATE_INTERVAL == 0 ) {
                for ( ReadsDownsampler<SAMRecord> perSampleDownsampler : perSampleDownsamplers.values() ) {
                    perSampleDownsampler.signalNoMoreReadsBefore(read);
                    processFinalizedAndPendingItems(perSampleDownsampler);
                }
            }

            prevRead = read;
        }

        if ( ! nestedSAMIterator.hasNext() ) {
            for ( ReadsDownsampler<SAMRecord> perSampleDownsampler : perSampleDownsamplers.values() ) {
                perSampleDownsampler.signalEndOfInput();
                if ( perSampleDownsampler.hasFinalizedItems() ) {
                    orderedDownsampledReadsCache.addAll(perSampleDownsampler.consumeFinalizedItems());
                }
            }
            earliestPendingRead = null;
            earliestPendingDownsampler = null;
        }

        return readyToReleaseReads();
    }

    private void updateEarliestPendingRead( ReadsDownsampler<SAMRecord> currentDownsampler ) {
        // If there is no recorded earliest pending read and this downsampler has pending items,
        // then this downsampler's first pending item becomes the new earliest pending read:
        if ( earliestPendingRead == null && currentDownsampler.hasPendingItems() ) {
            earliestPendingRead = currentDownsampler.peekPending();
            earliestPendingDownsampler = currentDownsampler;
        }
        // In all other cases, we only need to update the earliest pending read when the downsampler
        // associated with it experiences a change in its pending reads, since by assuming a sorted
        // read stream we're assured that each downsampler's earliest pending read will only increase
        // in genomic position over time.
        //
        // TODO: An occasional O(samples) linear search seems like a better option than keeping the downsamplers
        // TODO: sorted by earliest pending read, which would cost at least O(total_reads * (samples + log(samples))),
        // TODO: but need to verify this empirically.
        else if ( currentDownsampler == earliestPendingDownsampler &&
                  (! currentDownsampler.hasPendingItems() || readComparator.compare(currentDownsampler.peekPending(), earliestPendingRead) != 0) ) {

            earliestPendingRead = null;
            earliestPendingDownsampler = null;
            for ( ReadsDownsampler<SAMRecord> perSampleDownsampler : perSampleDownsamplers.values() ) {
                if ( perSampleDownsampler.hasPendingItems() &&
                     (earliestPendingRead == null || readComparator.compare(perSampleDownsampler.peekPending(), earliestPendingRead) < 0) ) {

                    earliestPendingRead = perSampleDownsampler.peekPending();
                    earliestPendingDownsampler = perSampleDownsampler;
                }
            }
        }
    }

    private void processFinalizedAndPendingItems( ReadsDownsampler<SAMRecord> currentDownsampler ) {
        if ( currentDownsampler.hasFinalizedItems() ) {
            orderedDownsampledReadsCache.addAll(currentDownsampler.consumeFinalizedItems());
        }
        updateEarliestPendingRead(currentDownsampler);
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    public void close() {
        nestedSAMIterator.close();
    }

    public Iterator<SAMRecord> iterator() {
        return this;
    }
}
