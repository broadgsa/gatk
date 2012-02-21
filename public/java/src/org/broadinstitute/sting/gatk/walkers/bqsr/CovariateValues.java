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
    private Comparable[] mismatches;
    private Comparable[] insertions;
    private Comparable[] deletions;

    public CovariateValues(Comparable[] mismatch, Comparable[] insertion, Comparable[] deletion) {
        this.mismatches = mismatch;
        this.insertions = insertion;
        this.deletions = deletion;
    }

    public Comparable[] getMismatches() {
        return mismatches;
    }

    public Comparable[] getInsertions() {
        return insertions;
    }

    public Comparable[] getDeletions() {
        return deletions;
    }

}
