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
    private final Long[][] mismatchesKeySet;
    private final Long[][] insertionsKeySet;
    private final Long[][] deletionsKeySet;

    private int nextCovariateIndex;

    public ReadCovariates(int readLength, int numberOfCovariates) {
        this.mismatchesKeySet = new Long[readLength][numberOfCovariates];
        this.insertionsKeySet = new Long[readLength][numberOfCovariates];
        this.deletionsKeySet = new Long[readLength][numberOfCovariates];
        this.nextCovariateIndex = 0;
    }

    public void addCovariate(CovariateValues covariate) {
        transposeCovariateValues(mismatchesKeySet, covariate.getMismatches());
        transposeCovariateValues(insertionsKeySet, covariate.getInsertions());
        transposeCovariateValues(deletionsKeySet, covariate.getDeletions());
        nextCovariateIndex++;
    }

    public Long[] getKeySet(final int readPosition, final EventType errorModel) {
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

    public Long[] getMismatchesKeySet(int readPosition) {
        return mismatchesKeySet[readPosition];
    }

    public Long[] getInsertionsKeySet(int readPosition) {
        return insertionsKeySet[readPosition];
    }

    public Long[] getDeletionsKeySet(int readPosition) {
        return deletionsKeySet[readPosition];
    }

    private void transposeCovariateValues(Long[][] keySet, Long[] covariateValues) {
        for (int i = 0; i < covariateValues.length; i++)
            keySet[i][nextCovariateIndex] = covariateValues[i];
    }

    /**
     * Testing routines
     */
    protected Long[][] getMismatchesKeySet() {
        return mismatchesKeySet;
    }

    protected Long[][] getInsertionsKeySet() {
        return insertionsKeySet;
    }

    protected Long[][] getDeletionsKeySet() {
        return deletionsKeySet;
    }
}
