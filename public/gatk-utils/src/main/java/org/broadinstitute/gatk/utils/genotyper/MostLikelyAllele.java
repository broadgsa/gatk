/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.genotyper;

import org.broadinstitute.gatk.utils.MathUtils;
import htsjdk.variant.variantcontext.Allele;

/**
 * Stores the most likely and second most likely alleles, along with a threshold
 * for assuming computing that a read is informative.
 *
 * If the difference between the most-likely allele and the next-most-likely allele is < INFORMATIVE_LIKELIHOOD_THRESHOLD
 * then the most likely allele is set to "no call", and isInformative will return false.  This constant can be
 * overridden simply by using one of the version of these calls that accepts informative threshold as an argument.
 *
 * For convenience, there are functions called getAlleleIfInformative that return either the most likely allele, or
 * NO_CALL if two or more alleles have likelihoods within INFORMATIVE_LIKELIHOOD_THRESHOLD of one another.
 *
 * By default empty allele maps will return NO_CALL, and allele maps with a single entry will return the
 * corresponding key
 *
 * User: depristo
 * Date: 3/24/13
 * Time: 1:39 PM
 */
public final class MostLikelyAllele {
    public static final double INFORMATIVE_LIKELIHOOD_THRESHOLD = 0.2;

    final Allele mostLikely;
    final Allele secondLikely;
    final double log10LikelihoodOfMostLikely;
    final double log10LikelihoodOfSecondBest;

    /**
     * Create a new MostLikelyAllele
     *
     * If there's a meaningful most likely allele, allele should be a real allele.  If none can be determined,
     * mostLikely should be a NO_CALL allele.
     *
     * @param mostLikely the most likely allele
     * @param secondMostLikely the most likely allele after mostLikely
     * @param log10LikelihoodOfMostLikely the log10 likelihood of the most likely allele
     * @param log10LikelihoodOfSecondBest the log10 likelihood of the next most likely allele (should be NEGATIVE_INFINITY if none is available)
     */
    public MostLikelyAllele(final Allele mostLikely, final Allele secondMostLikely, double log10LikelihoodOfMostLikely, double log10LikelihoodOfSecondBest) {
        if ( mostLikely == null ) throw new IllegalArgumentException("mostLikely allele cannot be null");
        if ( log10LikelihoodOfMostLikely != Double.NEGATIVE_INFINITY && ! MathUtils.goodLog10Probability(log10LikelihoodOfMostLikely) )
            throw new IllegalArgumentException("log10LikelihoodOfMostLikely must be either -Infinity or a good log10 prob but got " + log10LikelihoodOfMostLikely);
        if ( log10LikelihoodOfSecondBest != Double.NEGATIVE_INFINITY && ! MathUtils.goodLog10Probability(log10LikelihoodOfSecondBest) )
            throw new IllegalArgumentException("log10LikelihoodOfSecondBest must be either -Infinity or a good log10 prob but got " + log10LikelihoodOfSecondBest);
        if ( log10LikelihoodOfMostLikely < log10LikelihoodOfSecondBest )
            throw new IllegalArgumentException("log10LikelihoodOfMostLikely must be <= log10LikelihoodOfSecondBest but got " + log10LikelihoodOfMostLikely + " vs 2nd " + log10LikelihoodOfSecondBest);

        this.mostLikely = mostLikely;
        this.secondLikely = secondMostLikely;
        this.log10LikelihoodOfMostLikely = log10LikelihoodOfMostLikely;
        this.log10LikelihoodOfSecondBest = log10LikelihoodOfSecondBest;
    }

    public Allele getMostLikelyAllele() {
        return mostLikely;
    }

    public Allele getSecondMostLikelyAllele() {
        return secondLikely;
    }

    public double getLog10LikelihoodOfMostLikely() {
        return log10LikelihoodOfMostLikely;
    }

    public double getLog10LikelihoodOfSecondBest() {
        return log10LikelihoodOfSecondBest;
    }

    /**
     * @see #isInformative(double) with threshold of INFORMATIVE_LIKELIHOOD_THRESHOLD
     */
    public boolean isInformative() {
        return isInformative(INFORMATIVE_LIKELIHOOD_THRESHOLD);
    }

    /**
     * Was this allele selected from an object that was specifically informative about the allele?
     *
     * The calculation that implements this is whether the likelihood of the most likely allele is larger
     * than the second most likely by at least the log10ThresholdForInformative
     *
     * @return true if so, false if not
     */
    public boolean isInformative(final double log10ThresholdForInformative) {
        return getLog10LikelihoodOfMostLikely() - getLog10LikelihoodOfSecondBest() > log10ThresholdForInformative;
    }

    /**
     * @see #getAlleleIfInformative(double) with threshold of INFORMATIVE_LIKELIHOOD_THRESHOLD
     */
    public Allele getAlleleIfInformative() {
        return getAlleleIfInformative(INFORMATIVE_LIKELIHOOD_THRESHOLD);
    }

    /**
     * Get the most likely allele if isInformative(log10ThresholdForInformative) is true, or NO_CALL otherwise
     *
     * @param log10ThresholdForInformative a log10 threshold to determine if the most likely allele was informative
     * @return a non-null allele
     */
    public Allele getAlleleIfInformative(final double log10ThresholdForInformative) {
        return isInformative(log10ThresholdForInformative) ? getMostLikelyAllele() : Allele.NO_CALL;
    }
}
