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

package org.broadinstitute.gatk.engine.traversals;

import org.broadinstitute.gatk.utils.downsampling.Downsampler;
import org.broadinstitute.gatk.utils.downsampling.ReservoirDownsampler;
import org.broadinstitute.gatk.utils.sam.AlignmentStartComparator;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

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
    private final int maxCapacity;
    private ArrayList<GATKSAMRecord> undownsampledCache;
    private Downsampler<GATKSAMRecord> downsampler;

    private static final int UNDOWNSAMPLED_CACHE_MAX_INITIAL_SIZE = 10000;

    /**
     * Create a new empty ReadCache
     * @param maxCapacity the max capacity of the read cache.
     */
    public TAROrderedReadCache( final int maxCapacity ) {
        if ( maxCapacity < 0 ) throw new IllegalArgumentException("maxCapacity must be >= 0 but got " + maxCapacity);
        this.maxCapacity = maxCapacity;

        // The one we're not currently using will always be null:
        initializeUndownsampledCache();
        this.downsampler = null;
    }

    /**
     * Moves all reads over to the downsampler, causing it to be used from this point on. Should be called
     * when the undownsampledCache fills up and we need to start discarding reads. Since the
     * ReservoirDownsampler doesn't preserve relative ordering, pop operations become expensive
     * after this point, as they require a O(n log n) sort.
     */
    private void activateDownsampler() {
        downsampler = new ReservoirDownsampler<>(maxCapacity, false);
        downsampler.submit(undownsampledCache);
        undownsampledCache = null; // preferable to the O(n) clear() method
    }

    /**
     * Allocate the undownsampled cache used when we have fewer than maxCapacity items
     */
    private void initializeUndownsampledCache() {
        undownsampledCache = new ArrayList<>(Math.min(maxCapacity + 1, UNDOWNSAMPLED_CACHE_MAX_INITIAL_SIZE));
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
    public void add( final GATKSAMRecord read ) {
        if ( read == null ) throw new IllegalArgumentException("Read cannot be null");

        if ( downsampler != null ) {
            downsampler.submit(read);
        }
        else {
            undownsampledCache.add(read);

            // No more room in the undownsampledCache? Time to start downsampling
            if ( undownsampledCache.size() > maxCapacity ) {
                activateDownsampler();
            }
        }
    }

    /**
     * Add a collection of reads to this cache.  Assumed to be in sorted order w.r.t. the previously added reads and each other
     * @param reads a collection of reads to add
     */
    public void addAll( final List<GATKSAMRecord> reads ) {
        if ( reads == null ) throw new IllegalArgumentException("Reads cannot be null");
        for ( final GATKSAMRecord read : reads ) {
            add(read);
        }
    }

    /**
     * How many reads are currently in the cache?
     * @return a positive integer
     */
    public int size() {
        return downsampler != null ? downsampler.size() : undownsampledCache.size();
    }

    /**
     * How many reads were discarded since the last call to popCurrentReads
     *
     * @return number of items discarded during downsampling since last pop operation
     */
    public int getNumDiscarded() {
        return downsampler != null ? downsampler.getNumberOfDiscardedItems() : 0;
    }

    /**
     * Removes all reads currently in the cache, and returns them in sorted order (w.r.t. alignmentStart)
     *
     * Flushes this cache, so after this call the cache will contain no reads, and we'll be in the same
     * initial state as the constructor would put us in, with a non-null undownsampledCache and a null
     * downsampler.
     *
     * @return a list of GATKSAMRecords in this cache
     */
    public List<GATKSAMRecord> popCurrentReads() {
        final List<GATKSAMRecord> poppedReads;

        if ( downsampler == null ) {
            poppedReads = undownsampledCache;  // avoid making a copy here, since we're going to allocate a new cache
        }
        else {
            // If we triggered the downsampler, we need to sort the reads before returning them,
            // since the ReservoirDownsampler is not guaranteed to preserve relative ordering of items.
            // After consuming the downsampled items in this call to popCurrentReads(), we switch back
            // to using the undownsampledCache until we fill up again.
            poppedReads = downsampler.consumeFinalizedItems();  // avoid making a copy here
            Collections.sort(poppedReads, new AlignmentStartComparator());
            downsampler = null;
        }

        initializeUndownsampledCache();
        return poppedReads;
    }
}
