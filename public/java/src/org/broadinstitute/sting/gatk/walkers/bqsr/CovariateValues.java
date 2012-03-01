package org.broadinstitute.sting.gatk.walkers.bqsr;

/**
 * An object to hold the different covariate values for all bases in the read.
 *
 * Currently we have three different covariates for each read:
 *   - Mismatch
 *   - Insertion
 *   - Deletion
 *
 * @author Mauricio Carneiro
 * @since 2/8/12
 */
public class CovariateValues {
    private Object[] mismatches;
    private Object[] insertions;
    private Object[] deletions;

    public CovariateValues(Object[] mismatch, Object[] insertion, Object[] deletion) {
        this.mismatches = mismatch;
        this.insertions = insertion;
        this.deletions = deletion;
    }

    public Object[] getMismatches() {
        return mismatches;
    }

    public Object[] getInsertions() {
        return insertions;
    }

    public Object[] getDeletions() {
        return deletions;
    }

}
