/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.recalibration;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.LRUCache;
import org.broadinstitute.gatk.utils.recalibration.EventType;

/**
 * The object temporarily held by a read that describes all of it's covariates.
 *
 * In essence, this is an array of CovariateValues, but it also has some functionality to deal with the optimizations of the NestedHashMap
 *
 * @author Mauricio Carneiro
 * @since 2/8/12
 */
public class ReadCovariates {
    private final static Logger logger = Logger.getLogger(ReadCovariates.class);

    /**
     * How big should we let the LRU cache grow
     */
    private static final int LRU_CACHE_SIZE = 500;

    /**
     * Use an LRU cache to keep cache of keys (int[][][]) arrays for each read length we've seen.
     * The cache allows us to avoid the expense of recreating these arrays for every read.  The LRU
     * keeps the total number of cached arrays to less than LRU_CACHE_SIZE.
     *
     * This is a thread local variable, so the total memory required may grow to N_THREADS x LRU_CACHE_SIZE
     */
    private final static ThreadLocal<LRUCache<Integer, int[][][]>> keysCache = new ThreadLocal<LRUCache<Integer, int[][][]>>() {
        @Override protected LRUCache<Integer, int[][][]> initialValue() {
            return new LRUCache<Integer, int[][][]>(LRU_CACHE_SIZE);
        }
    };

    /**
     * The keys cache is only valid for a single covariate count.  Normally this will remain constant for the analysis.
     * If running multiple analyses (or the unit test suite), it's necessary to clear the cache.
     */
    public static void clearKeysCache() {
        keysCache.remove();
    }

    /**
     * Our keys, indexed by event type x read length x covariate
     */
    private final int[][][] keys;

    /**
     * The index of the current covariate, used by addCovariate
     */
    private int currentCovariateIndex = 0;

    public ReadCovariates(final int readLength, final int numberOfCovariates) {
        final LRUCache<Integer, int[][][]> cache = keysCache.get();
        final int[][][] cachedKeys = cache.get(readLength);
        if ( cachedKeys == null ) {
            // There's no cached value for read length so we need to create a new int[][][] array
            if ( logger.isDebugEnabled() ) logger.debug("Keys cache miss for length " + readLength + " cache size " + cache.size());
            keys = new int[EventType.values().length][readLength][numberOfCovariates];
            cache.put(readLength, keys);
        } else {
            keys = cachedKeys;
        }
    }

    public void setCovariateIndex(final int index) {
        currentCovariateIndex = index;
    }

    /**
     * Update the keys for mismatch, insertion, and deletion for the current covariate at read offset
     *
     * NOTE: no checks are performed on the number of covariates, for performance reasons.  If the count increases
     * after the keysCache has been accessed, this method will throw an ArrayIndexOutOfBoundsException.  This currently
     * only occurs in the testing harness, and we don't anticipate that it will become a part of normal runs.
     *
     * @param mismatch the mismatch key value
     * @param insertion the insertion key value
     * @param deletion the deletion key value
     * @param readOffset the read offset, must be >= 0 and <= the read length used to create this ReadCovariates
     */
    public void addCovariate(final int mismatch, final int insertion, final int deletion, final int readOffset) {
        keys[EventType.BASE_SUBSTITUTION.ordinal()][readOffset][currentCovariateIndex] = mismatch;
        keys[EventType.BASE_INSERTION.ordinal()][readOffset][currentCovariateIndex] = insertion;
        keys[EventType.BASE_DELETION.ordinal()][readOffset][currentCovariateIndex] = deletion;
    }

    /**
     * Get the keys for all covariates at read position for error model
     *
     * @param readPosition
     * @param errorModel
     * @return
     */
    public int[] getKeySet(final int readPosition, final EventType errorModel) {
        return keys[errorModel.ordinal()][readPosition];
    }

    public int[][] getKeySet(final EventType errorModel) {
        return keys[errorModel.ordinal()];
    }

    // ----------------------------------------------------------------------
    //
    // routines for testing
    //
    // ----------------------------------------------------------------------

    protected int[][] getMismatchesKeySet() { return getKeySet(EventType.BASE_SUBSTITUTION); }
    protected int[][] getInsertionsKeySet() { return getKeySet(EventType.BASE_INSERTION); }
    protected int[][] getDeletionsKeySet() { return getKeySet(EventType.BASE_DELETION); }

    protected int[] getMismatchesKeySet(final int readPosition) {
        return getKeySet(readPosition, EventType.BASE_SUBSTITUTION);
    }

    protected int[] getInsertionsKeySet(final int readPosition) {
        return getKeySet(readPosition, EventType.BASE_INSERTION);
    }

    protected int[] getDeletionsKeySet(final int readPosition) {
        return getKeySet(readPosition, EventType.BASE_DELETION);
    }
}
