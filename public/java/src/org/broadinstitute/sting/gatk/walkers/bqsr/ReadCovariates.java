package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.BitSet;

/**
 * The object temporarily held by a read that describes all of it's covariates.
 *
 * In essence, this is an array of CovariateValues, but it also has some functionality to deal with the optimizations of the NestedHashMap
 *
 * @author Mauricio Carneiro
 * @since 2/8/12
 */
public class ReadCovariates {
    private final BitSet[][] mismatchesKeySet;
    private final BitSet[][] insertionsKeySet;
    private final BitSet[][] deletionsKeySet;

    private int nextCovariateIndex;

    public ReadCovariates(int readLength, int numberOfCovariates) {
        this.mismatchesKeySet = new BitSet[readLength][numberOfCovariates];
        this.insertionsKeySet = new BitSet[readLength][numberOfCovariates];
        this.deletionsKeySet = new BitSet[readLength][numberOfCovariates];
        this.nextCovariateIndex = 0;
    }

    public void addCovariate(CovariateValues covariate) {
        transposeCovariateValues(mismatchesKeySet, covariate.getMismatches());
        transposeCovariateValues(insertionsKeySet, covariate.getInsertions());
        transposeCovariateValues(deletionsKeySet, covariate.getDeletions());
        nextCovariateIndex++;
    }

    public BitSet[] getKeySet(final int readPosition, final EventType errorModel) {
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

    public BitSet[] getMismatchesKeySet(int readPosition) {
        return mismatchesKeySet[readPosition];
    }

    public BitSet[] getInsertionsKeySet(int readPosition) {
        return insertionsKeySet[readPosition];
    }

    public BitSet[] getDeletionsKeySet(int readPosition) {
        return deletionsKeySet[readPosition];
    }

    private void transposeCovariateValues(BitSet[][] keySet, BitSet[] covariateValues) {
        for (int i = 0; i < covariateValues.length; i++)
            keySet[i][nextCovariateIndex] = covariateValues[i];
    }
}
