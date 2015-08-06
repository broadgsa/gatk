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

package org.broadinstitute.gatk.engine.iterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;

/**
 * Iterates through a list of elements, tracking the number of elements it has seen.
 * @author hanna
 * @version 0.1
 */
public class PositionTrackingIterator implements GATKSAMIterator {
    /**
     * The iterator being tracked.
     */
    private CloseableIterator<SAMRecord> iterator;

    /**
     * Current position within the tracked iterator.
     */
    private long position;

    /**
     * Retrieves the current position of the iterator.  The 'current position' of the iterator is defined as
     * the coordinate of the read that will be returned if next() is called.
     * @return The current position of the iterator.
     */
    public long getPosition() {
        return position;
    }

    /**
     * Create a new iterator wrapping the given position, assuming that the reader is <code>position</code> reads
     * into the sequence.
     * @param iterator Iterator to wraps.
     * @param position Non-negative position where the iterator currently sits.
     */
    public PositionTrackingIterator(CloseableIterator<SAMRecord> iterator, long position ) {
        this.iterator = iterator;
        this.position = position;
    }

    /**
     * {@inheritDoc}
     */
    public boolean hasNext() {
        return iterator.hasNext();
    }

    /**
     * Try to get the next read in the list.  If a next read is available, increment the position.
     * @return next read in the list, if available.
     */
    public SAMRecord next() {
        try {
            return iterator.next();
        }
        finally {
            position++;
        }
    }

    /**
     * {@inheritDoc}
     */
    public GATKSAMIterator iterator() {
        return this;
    }

    /**
     * {@inheritDoc}
     */
    public void close() {
        iterator.close();
    }

    /**
     * {@inheritDoc}
     */
    public void remove() { throw new UnsupportedOperationException("Cannot remove from a GATKSAMIterator"); }
}
