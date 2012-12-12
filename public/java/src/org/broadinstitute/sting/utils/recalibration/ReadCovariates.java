package org.broadinstitute.sting.utils.recalibration;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.LRUCache;

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
            logger.info("Keys cache miss for length " + readLength + " cache size " + cache.size());
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
     * @param mismatch the mismatch key value
     * @param insertion the insertion key value
     * @param deletion the deletion key value
     * @param readOffset the read offset, must be >= 0 and <= the read length used to create this ReadCovariates
     */
    public void addCovariate(final int mismatch, final int insertion, final int deletion, final int readOffset) {
        keys[EventType.BASE_SUBSTITUTION.index][readOffset][currentCovariateIndex] = mismatch;
        keys[EventType.BASE_INSERTION.index][readOffset][currentCovariateIndex] = insertion;
        keys[EventType.BASE_DELETION.index][readOffset][currentCovariateIndex] = deletion;
    }

    /**
     * Get the keys for all covariates at read position for error model
     *
     * @param readPosition
     * @param errorModel
     * @return
     */
    public int[] getKeySet(final int readPosition, final EventType errorModel) {
        return keys[errorModel.index][readPosition];
    }

    public int[][] getKeySet(final EventType errorModel) {
        return keys[errorModel.index];
    }

    // ----------------------------------------------------------------------
    //
    // routines for testing
    //
    // ----------------------------------------------------------------------

    protected int[][] getMismatchesKeySet() {
        return keys[EventType.BASE_SUBSTITUTION.index];
    }

    protected int[][] getInsertionsKeySet() {
        return keys[EventType.BASE_INSERTION.index];
    }

    protected int[][] getDeletionsKeySet() {
        return keys[EventType.BASE_DELETION.index];
    }

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
