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
    private final Long[] mismatches;
    private final Long[] insertions;
    private final Long[] deletions;

    public CovariateValues(Long[] mismatch, Long[] insertion, Long[] deletion) {
        this.mismatches = mismatch;
        this.insertions = insertion;
        this.deletions = deletion;
    }

    public Long[] getMismatches() {
        return mismatches;
    }

    public Long[] getInsertions() {
        return insertions;
    }

    public Long[] getDeletions() {
        return deletions;
    }

}
