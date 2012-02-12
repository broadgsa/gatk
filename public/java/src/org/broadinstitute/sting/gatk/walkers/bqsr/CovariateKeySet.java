package org.broadinstitute.sting.gatk.walkers.bqsr;

/**
 * The object temporarily held by a read that describes all of it's covariates. 
 * 
 * In essence, this is an array of CovariateValues, but it also has some functionality to deal with the optimizations of the NestedHashMap
 *
 * @author Mauricio Carneiro
 * @since 2/8/12
 */
public class CovariateKeySet {
    private Object[][] mismatchesKeySet;
    private Object[][] insertionsKeySet;
    private Object[][]  deletionsKeySet;

    private int nextCovariateIndex;
    
    public final static String mismatchesCovariateName = "M";
    public final static String insertionsCovariateName = "I";
    public final static String  deletionsCovariateName = "D";

    public CovariateKeySet(int readLength, int numberOfCovariates) {
        numberOfCovariates++;                                               // +1 because we are adding the mismatch covariate (to comply with the molten table format)
        this.mismatchesKeySet = new Object[readLength][numberOfCovariates]; 
        this.insertionsKeySet = new Object[readLength][numberOfCovariates];
        this.deletionsKeySet  = new Object[readLength][numberOfCovariates];
        initializeCovariateKeySet(this.mismatchesKeySet, mismatchesCovariateName);
        initializeCovariateKeySet(this.insertionsKeySet, insertionsCovariateName);
        initializeCovariateKeySet(this.deletionsKeySet,  deletionsCovariateName);
        this.nextCovariateIndex = 0;
    }
    
    public void addCovariate(CovariateValues covariate) {
        transposeCovariateValues(mismatchesKeySet, covariate.getMismatches());
        transposeCovariateValues(insertionsKeySet, covariate.getInsertions());
        transposeCovariateValues(deletionsKeySet,  covariate.getDeletions());
        nextCovariateIndex++;
    }
    
    public Object[] getMismatchesKeySet(int readPosition) {
        return mismatchesKeySet[readPosition];
    }

    public Object[] getInsertionsKeySet(int readPosition) {
        return insertionsKeySet[readPosition];
    }

    public Object[] getDeletionsKeySet(int readPosition) {
        return deletionsKeySet[readPosition];
    }

    private void transposeCovariateValues (Object [][] keySet, Object [] covariateValues) {
        for (int i=0; i<covariateValues.length; i++) 
            keySet[i][nextCovariateIndex] = covariateValues[i];        
    }
    
    private void initializeCovariateKeySet (Object[][] keySet, String covariateName) {
        int readLength = keySet.length;
        int lastCovariateIndex = keySet[0].length - 1;
        for (int i = 0; i < readLength; i++) 
            keySet[i][lastCovariateIndex] = covariateName;
    }
}
