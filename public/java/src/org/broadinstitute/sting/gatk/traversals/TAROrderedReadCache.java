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

package org.broadinstitute.sting.gatk.traversals;

import org.broadinstitute.sting.gatk.downsampling.Downsampler;
import org.broadinstitute.sting.gatk.downsampling.ReservoirDownsampler;
import org.broadinstitute.sting.utils.sam.AlignmentStartComparator;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Subsystem to track a list of all reads currently live in the TraverseActiveRegions system,
 * while limiting the total number of reads to a maximum capacity.
 *
 * User: depristo
 * Date: 4/7/13
 * Time: 11:23 AM
 */
public class TAROrderedReadCache {
    final int maxCapacity;
    final Downsampler<GATKSAMRecord> downsampler;

    /**
     * Create a new empty ReadCache
     * @param maxCapacity the max capacity of the read cache.
     */
    public TAROrderedReadCache(int maxCapacity) {
        if ( maxCapacity < 0 ) throw new IllegalArgumentException("maxCapacity must be >= 0 but got " + maxCapacity);
        this.maxCapacity = maxCapacity;
        this.downsampler = new ReservoirDownsampler<GATKSAMRecord>(maxCapacity);
    }

    /**
     * What's the maximum number of reads we'll store in the cache?
     * @return a positive integer
     */
    public int getMaxCapacity() {
        return maxCapacity;
    }

    /**
     * Add a single read to this cache.  Assumed to be in sorted order w.r.t. the previously added reads
     * @param read a read to add
     */
    public void add(final GATKSAMRecord read) {
        if ( read == null ) throw new IllegalArgumentException("Read cannot be null");
        downsampler.submit(read);
    }

    /**
     * Add a collection of reads to this cache.  Assumed to be in sorted order w.r.t. the previously added reads and each other
     * @param reads a collection of reads to add
     */
    public void addAll(final List<GATKSAMRecord> reads) {
        if ( reads == null ) throw new IllegalArgumentException("Reads cannot be null");
        downsampler.submit(reads);
    }

    /**
     * How many reads are currently in the cache?
     * @return a positive integer
     */
    public int size() {
        return downsampler.size();
    }

    /**
     * How many reads were discarded since the last call to popCurrentReads
     * @return
     */
    public int getNumDiscarded() {
        return downsampler.getNumberOfDiscardedItems();
    }

    /**
     * Removes all reads currently in the cache, and returns them in sorted order (w.r.t. alignmentStart)
     *
     * Flushes this cache, so after this call the cache will contain no reads and all downsampling stats will
     * be reset.
     *
     * @return a list of GATKSAMRecords in this cache
     */
    public List<GATKSAMRecord> popCurrentReads() {
        final List<GATKSAMRecord> maybeUnordered = downsampler.consumeFinalizedItems();

        final List<GATKSAMRecord> ordered;
        if ( downsampler.getNumberOfDiscardedItems() == 0 ) {
            // haven't discarded anything, so the reads are ordered properly
            ordered = maybeUnordered;
        } else {
            // we need to sort these damn things: O(n log n)
            ordered = new ArrayList<GATKSAMRecord>(maybeUnordered);
            Collections.sort(ordered, new AlignmentStartComparator());
        }

        // reset the downsampler stats so getNumberOfDiscardedItems is 0
        downsampler.reset();
        return ordered;
    }
}
