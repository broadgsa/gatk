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
    private final long[] mismatches;
    private final long[] insertions;
    private final long[] deletions;

    public CovariateValues(final long[] mismatch, final long[] insertion, final long[] deletion) {
        this.mismatches = mismatch;
        this.insertions = insertion;
        this.deletions = deletion;
    }

    public long[] getMismatches() {
        return mismatches;
    }

    public long[] getInsertions() {
        return insertions;
    }

    public long[] getDeletions() {
        return deletions;
    }

}
