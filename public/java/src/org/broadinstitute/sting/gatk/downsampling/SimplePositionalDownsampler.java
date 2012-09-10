/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.downsampling;

import net.sf.samtools.SAMRecord;

import java.util.*;

/**
 * Simple Positional Downsampler: Downsample each stack of reads at each alignment start to a size <= a target coverage
 * using a Reservoir downsampler. Stores only O(target coverage) reads in memory at any given time.
 *
 * @author David Roazen
 */
public class SimplePositionalDownsampler<T extends SAMRecord> implements ReadsDownsampler<T> {

    private int targetCoverage;

    private ReservoirDownsampler<T> reservoir;

    private int currentContigIndex;

    private int currentAlignmentStart;

    private boolean positionEstablished;

    private boolean unmappedReadsReached;

    private ArrayList<T> finalizedReads;

    private int numDiscardedItems;

    /**
     * Construct a SimplePositionalDownsampler
     *
     * @param targetCoverage Maximum number of reads that may share any given alignment start position
     */
    public SimplePositionalDownsampler( int targetCoverage ) {
        this.targetCoverage = targetCoverage;
        reservoir = new ReservoirDownsampler<T>(targetCoverage);
        finalizedReads = new ArrayList<T>();
        clear();
        reset();
    }

    public void submit( T newRead ) {
        updatePositionalState(newRead);

        if ( unmappedReadsReached ) {    // don't downsample the unmapped reads at the end of the stream
            finalizedReads.add(newRead);
        }
        else {
            int reservoirPreviouslyDiscardedItems = reservoir.getNumberOfDiscardedItems();
            reservoir.submit(newRead);
            numDiscardedItems += reservoir.getNumberOfDiscardedItems() - reservoirPreviouslyDiscardedItems;
        }
    }

    public void submit( Collection<T> newReads ) {
        for ( T read : newReads ) {
            submit(read);
        }
    }

    public boolean hasFinalizedItems() {
        return finalizedReads.size() > 0;
    }

    public List<T> consumeFinalizedItems() {
        // pass by reference rather than make a copy, for speed
        List<T> toReturn = finalizedReads;
        finalizedReads = new ArrayList<T>();
        return toReturn;
    }

    public boolean hasPendingItems() {
        return reservoir.hasFinalizedItems();
    }

    public T peekFinalized() {
        return finalizedReads.isEmpty() ? null : finalizedReads.get(0);
    }

    public T peekPending() {
        return reservoir.peekFinalized();
    }

    public int getNumberOfDiscardedItems() {
        return numDiscardedItems;
    }

    public void signalEndOfInput() {
        finalizeReservoir();
    }

    public void clear() {
        reservoir.clear();
        reservoir.reset();
        finalizedReads.clear();
        positionEstablished = false;
        unmappedReadsReached = false;
    }

    public void reset() {
        numDiscardedItems = 0;
    }

    public boolean requiresCoordinateSortOrder() {
        return true;
    }

    public void signalNoMoreReadsBefore( T read ) {
        updatePositionalState(read);
    }

    private void updatePositionalState( T newRead ) {
        if ( readIsPastCurrentPosition(newRead) ) {
            if ( reservoir.hasFinalizedItems() ) {
                finalizeReservoir();
            }

            setCurrentPosition(newRead);

            if ( newRead.getReadUnmappedFlag() ) {
                unmappedReadsReached = true;
            }
        }
    }

    private void setCurrentPosition( T read ) {
        currentContigIndex = read.getReferenceIndex();
        currentAlignmentStart = read.getAlignmentStart();
        positionEstablished = true;
    }

    private boolean readIsPastCurrentPosition( T read ) {
        return ! positionEstablished ||
               read.getReferenceIndex() > currentContigIndex ||
               read.getAlignmentStart() > currentAlignmentStart ||
               (read.getReadUnmappedFlag() && ! unmappedReadsReached);
    }

    private void finalizeReservoir() {
        finalizedReads.addAll(reservoir.consumeFinalizedItems());
        reservoir.reset();
    }
}
