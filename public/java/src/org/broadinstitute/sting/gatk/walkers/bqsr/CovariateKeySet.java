package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.BitSet;
import java.util.HashMap;

/**
 * The object temporarily held by a read that describes all of it's covariates.
 *
 * In essence, this is an array of CovariateValues, but it also has some functionality to deal with the optimizations of the NestedHashMap
 *
 * @author Mauricio Carneiro
 * @since 2/8/12
 */
public class CovariateKeySet {
    private BitSet[][] mismatchesKeySet;
    private BitSet[][] insertionsKeySet;
    private BitSet[][] deletionsKeySet;

    private int nextCovariateIndex;

    //    private static String mismatchesCovariateName = "M";
    //    private static String insertionsCovariateName = "I";
    //    private static String deletionsCovariateName  = "D";
    //
    //    private static BitSet mismatchesCovariateBitSet = BitSetUtils.bitSetFrom(0);
    //    private static BitSet insertionsCovariateBitSet = BitSetUtils.bitSetFrom(1);
    //    private static BitSet deletionsCovariateBitSet = BitSetUtils.bitSetFrom(2);

    private static HashMap<String, RecalDataManager.BaseRecalibrationType> nameToType = new HashMap<String, RecalDataManager.BaseRecalibrationType>();
    private static HashMap<BitSet, String> bitSetToName = new HashMap<BitSet, String>();

    public CovariateKeySet(int readLength, int numberOfCovariates) {
        //        numberOfCovariates++;                                               // +1 because we are adding the mismatch covariate (to comply with the molten table format)
        this.mismatchesKeySet = new BitSet[readLength][numberOfCovariates];
        this.insertionsKeySet = new BitSet[readLength][numberOfCovariates];
        this.deletionsKeySet = new BitSet[readLength][numberOfCovariates];
        //        initializeCovariateKeySet(this.mismatchesKeySet, mismatchesCovariateBitSet);
        //        initializeCovariateKeySet(this.insertionsKeySet, insertionsCovariateBitSet);
        //        initializeCovariateKeySet(this.deletionsKeySet, deletionsCovariateBitSet);
        this.nextCovariateIndex = 0;

        //        nameToType.put(mismatchesCovariateName, RecalDataManager.BaseRecalibrationType.BASE_SUBSTITUTION);
        //        nameToType.put(insertionsCovariateName, RecalDataManager.BaseRecalibrationType.BASE_INSERTION);
        //        nameToType.put(deletionsCovariateName,  RecalDataManager.BaseRecalibrationType.BASE_DELETION);
        //
        //        bitSetToName.put(BitSetUtils.bitSetFrom(0), mismatchesCovariateName);
        //        bitSetToName.put(BitSetUtils.bitSetFrom(1), insertionsCovariateName);
        //        bitSetToName.put(BitSetUtils.bitSetFrom(2), deletionsCovariateName);
    }

    public void addCovariate(CovariateValues covariate) {
        transposeCovariateValues(mismatchesKeySet, covariate.getMismatches());
        transposeCovariateValues(insertionsKeySet, covariate.getInsertions());
        transposeCovariateValues(deletionsKeySet, covariate.getDeletions());
        nextCovariateIndex++;
    }

    public static RecalDataManager.BaseRecalibrationType errorModelFrom(final String modelString) {
        if (!nameToType.containsKey(modelString))
            throw new ReviewedStingException("Unrecognized Base Recalibration model string: " + modelString);
        return nameToType.get(modelString);
    }

    public static String eventNameFrom(final BitSet bitSet) {
        if (!bitSetToName.containsKey(bitSet))
            throw new ReviewedStingException("Unrecognized Event Type BitSet: " + bitSet);
        return bitSetToName.get(bitSet);
    }

    public BitSet[] getKeySet(final int readPosition, final RecalDataManager.BaseRecalibrationType errorModel) {
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

    private void initializeCovariateKeySet(BitSet[][] keySet, BitSet covariateName) {
        int readLength = keySet.length;
        int lastCovariateIndex = keySet[0].length - 1;
        for (int i = 0; i < readLength; i++)
            keySet[i][lastCovariateIndex] = covariateName;
    }
}
