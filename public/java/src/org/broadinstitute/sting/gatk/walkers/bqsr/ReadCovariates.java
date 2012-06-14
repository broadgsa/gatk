package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

/**
 * The object temporarily held by a read that describes all of it's covariates.
 *
 * In essence, this is an array of CovariateValues, but it also has some functionality to deal with the optimizations of the NestedHashMap
 *
 * @author Mauricio Carneiro
 * @since 2/8/12
 */
public class ReadCovariates {
    private final long[][] mismatchesKeySet;
    private final long[][] insertionsKeySet;
    private final long[][] deletionsKeySet;

    private int nextCovariateIndex;

    public ReadCovariates(int readLength, int numberOfCovariates) {
        this.mismatchesKeySet = new long[readLength][numberOfCovariates];
        this.insertionsKeySet = new long[readLength][numberOfCovariates];
        this.deletionsKeySet = new long[readLength][numberOfCovariates];
        this.nextCovariateIndex = 0;
    }

    public void reset() {
        nextCovariateIndex = 0;
    }

    public void addCovariate(CovariateValues covariate) {
        transposeCovariateValues(mismatchesKeySet, covariate.getMismatches());
        transposeCovariateValues(insertionsKeySet, covariate.getInsertions());
        transposeCovariateValues(deletionsKeySet, covariate.getDeletions());
        nextCovariateIndex++;
    }

    public long[] getKeySet(final int readPosition, final EventType errorModel) {
        switch (errorModel) {
            case BASE_SUBSTITUTION:
                return getMismatchesKeySet(readPosition);
            case BASE_INSERTION:
                return getInsertionsKeySet(readPosition);
            case BASE_DELETION:
                return getDeletionsKeySet(readPosition);
            default:
                throw new ReviewedStingException("Unrecognized Base Recalibration type: " + errorModel);
        }
    }

    public long[] getMismatchesKeySet(final int readPosition) {
        return mismatchesKeySet[readPosition];
    }

    public long[] getInsertionsKeySet(final int readPosition) {
        return insertionsKeySet[readPosition];
    }

    public long[] getDeletionsKeySet(final int readPosition) {
        return deletionsKeySet[readPosition];
    }

    private void transposeCovariateValues(final long[][] keySet, final long[] covariateValues) {
        for (int i = 0; i < covariateValues.length; i++)
            keySet[i][nextCovariateIndex] = covariateValues[i];
    }

    /**
     * Testing routines
     */
    protected long[][] getMismatchesKeySet() {
        return mismatchesKeySet;
    }

    protected long[][] getInsertionsKeySet() {
        return insertionsKeySet;
    }

    protected long[][] getDeletionsKeySet() {
        return deletionsKeySet;
    }
}
