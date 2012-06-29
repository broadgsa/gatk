package org.broadinstitute.sting.gatk.walkers.bqsr;

/**
 * The object temporarily held by a read that describes all of it's covariates.
 *
 * In essence, this is an array of CovariateValues, but it also has some functionality to deal with the optimizations of the NestedHashMap
 *
 * @author Mauricio Carneiro
 * @since 2/8/12
 */
public class ReadCovariates {
    private final int[][][] keys;

    private int currentCovariateIndex = 0;

    public ReadCovariates(final int readLength, final int numberOfCovariates) {
        keys = new int[EventType.values().length][readLength][numberOfCovariates];
    }

    public void setCovariateIndex(final int index) {
        currentCovariateIndex = index;
    }

    public void addCovariate(final int mismatch, final int insertion, final int deletion, final int readOffset) {
        keys[EventType.BASE_SUBSTITUTION.index][readOffset][currentCovariateIndex] = mismatch;
        keys[EventType.BASE_INSERTION.index][readOffset][currentCovariateIndex] = insertion;
        keys[EventType.BASE_DELETION.index][readOffset][currentCovariateIndex] = deletion;
    }

    public int[] getKeySet(final int readPosition, final EventType errorModel) {
        return keys[errorModel.index][readPosition];
    }

    public int[][] getKeySet(final EventType errorModel) {
        return keys[errorModel.index];
    }

    public int[] getMismatchesKeySet(final int readPosition) {
        return keys[EventType.BASE_SUBSTITUTION.index][readPosition];
    }

    public int[] getInsertionsKeySet(final int readPosition) {
        return keys[EventType.BASE_INSERTION.index][readPosition];
    }

    public int[] getDeletionsKeySet(final int readPosition) {
        return keys[EventType.BASE_DELETION.index][readPosition];
    }

    /**
     * Testing routines
     */
    protected int[][] getMismatchesKeySet() {
        return keys[EventType.BASE_SUBSTITUTION.index];
    }

    protected int[][] getInsertionsKeySet() {
        return keys[EventType.BASE_INSERTION.index];
    }

    protected int[][] getDeletionsKeySet() {
        return keys[EventType.BASE_DELETION.index];
    }
}
