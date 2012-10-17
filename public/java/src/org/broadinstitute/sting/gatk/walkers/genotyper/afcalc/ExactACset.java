package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.broadinstitute.sting.utils.MathUtils;

import java.util.Arrays;

/**
* Created with IntelliJ IDEA.
* User: depristo
* Date: 10/5/12
* Time: 2:53 PM
* To change this template use File | Settings | File Templates.
*/ // This class represents a column in the Exact AC calculation matrix
public final class ExactACset {
    // the counts of the various alternate alleles which this column represents
    private final ExactACcounts ACcounts;

    // the column of the matrix
    private final double[] log10Likelihoods;

    int sum = -1;

    public ExactACset(final int size, final ExactACcounts ACcounts) {
        this.ACcounts = ACcounts;
        log10Likelihoods = new double[size];
        Arrays.fill(log10Likelihoods, Double.NEGATIVE_INFINITY);
    }

    /**
     * sum of all the non-reference alleles
     */
    public int getACsum() {
        if ( sum == -1 )
            sum = (int)MathUtils.sum(getACcounts().getCounts());
        return sum;
    }

    public boolean equals(Object obj) {
        return (obj instanceof ExactACset) && getACcounts().equals(((ExactACset)obj).getACcounts());
    }

    public ExactACcounts getACcounts() {
        return ACcounts;
    }

    public double[] getLog10Likelihoods() {
        return log10Likelihoods;
    }
}
