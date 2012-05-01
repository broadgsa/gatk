package org.broadinstitute.sting.gatk.walkers.bqsr;

import java.util.BitSet;

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
    private final BitSet[] mismatches;
    private final BitSet[] insertions;
    private final BitSet[] deletions;

    public CovariateValues(BitSet[] mismatch, BitSet[] insertion, BitSet[] deletion) {
        this.mismatches = mismatch;
        this.insertions = insertion;
        this.deletions = deletion;
    }

    public BitSet[] getMismatches() {
        return mismatches;
    }

    public BitSet[] getInsertions() {
        return insertions;
    }

    public BitSet[] getDeletions() {
        return deletions;
    }

}
