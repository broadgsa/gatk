/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import com.google.java.contract.Ensures;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Dec 14, 2011
 *
 * Useful helper class to communicate the results of the allele frequency calculation
 *
 * TODO -- WHAT IS THE CONTRACT ON MAP AC AND P NON REF?
 */
class AFCalcResultTracker {
    protected static final double VALUE_NOT_CALCULATED = Double.NEGATIVE_INFINITY;

    // These variables are intended to contain the MLE and MAP (and their corresponding allele counts) of the site over all alternate alleles
    protected double log10MLE;
    protected double log10MAP;
    private final int[] alleleCountsOfMLE;
    private final int[] alleleCountsOfMAP;

    // The posteriors seen, not including that of AF=0
    private static final int LIKELIHOODS_CACHE_SIZE = 5000;
    private final double[] log10LikelihoodsMatrixValues = new double[LIKELIHOODS_CACHE_SIZE];
    private int currentLikelihoodsCacheIndex = 0;
    protected Double log10LikelihoodsMatrixSum = null;

    // These variables are intended to contain the likelihood/posterior probability for the site's being monomorphic (i.e. AF=0 for all alternate alleles)
    private double log10LikelihoodOfAFzero;
    private double log10PosteriorOfAFzero;
    private int[] AClimits;

    int nEvaluations = 0;

    /**
     * The list of alleles actually used in computing the AF
     */
    private List<Allele> allelesUsedInGenotyping = null;

    /**
     * Create a results object capability of storing results for calls with up to maxAltAlleles
     *
     * @param maxAltAlleles an integer >= 1
     */
    public AFCalcResultTracker(final int maxAltAlleles) {
        if ( maxAltAlleles < 1 ) throw new IllegalArgumentException("maxAltAlleles must be >= 0, saw " + maxAltAlleles);

        alleleCountsOfMLE = new int[maxAltAlleles];
        alleleCountsOfMAP = new int[maxAltAlleles];

        reset();
    }

    /**
     * Returns a vector with maxAltAlleles values containing AC values at the MLE
     *
     * The values of the ACs for this call are stored in the getAllelesUsedInGenotyping order,
     * starting from index 0 (i.e., the first alt allele is at 0).  The vector is always
     * maxAltAlleles in length, and so only the first getAllelesUsedInGenotyping.size() - 1 values
     * are meaningful.
     *
     * @return a vector with allele counts, not all of which may be meaningful
     */
    @Ensures("result != null")
    public int[] getAlleleCountsOfMLE() {
        return alleleCountsOfMLE;
    }

    /**
     * Returns a vector with maxAltAlleles values containing AC values at the MAP
     *
     * @see #getAlleleCountsOfMLE() for the encoding of results in this vector
     *
     * @return a non-null vector of ints
     */
    @Ensures("result != null")
    public int[] getAlleleCountsOfMAP() {
        return alleleCountsOfMAP;
    }

    /**
     * Returns the likelihoods summed across all AC values for AC > 0
     *
     * @return
     */
    public double getLog10LikelihoodOfAFNotZero() {
        if ( log10LikelihoodsMatrixSum == null ) {
            if ( currentLikelihoodsCacheIndex == 0 ) // there's nothing to sum up, so make the sum equal to the smallest thing we have
                log10LikelihoodsMatrixSum = MathUtils.LOG10_P_OF_ZERO;
            else
                log10LikelihoodsMatrixSum = MathUtils.log10sumLog10(log10LikelihoodsMatrixValues, 0, currentLikelihoodsCacheIndex);
        }
        return log10LikelihoodsMatrixSum;
    }

    public double getLog10LikelihoodOfAFNotZero(final boolean capAt0) {
        return Math.min(getLog10LikelihoodOfAFNotZero(), capAt0 ? 0.0 : Double.POSITIVE_INFINITY);
    }

    /**
     * TODO -- eric what is this supposed to return?  my unit tests don't do what I think they should
     *
     * @return
     */
    public double getLog10LikelihoodOfAFzero() {
        return log10LikelihoodOfAFzero;
    }

    /**
     * TODO -- eric what is this supposed to return?  my unit tests don't do what I think they should
     *
     * @return
     */
    public double getLog10PosteriorOfAFzero() {
        return log10PosteriorOfAFzero;
    }

    protected AFCalcResult toAFCalcResult(final double[] log10PriorsByAC) {
        final int [] subACOfMLE = Arrays.copyOf(alleleCountsOfMLE, allelesUsedInGenotyping.size() - 1);
        final double[] log10Likelihoods = new double[]{getLog10LikelihoodOfAFzero(), getLog10LikelihoodOfAFNotZero(true)};
        final double[] log10Priors = MathUtils.normalizeFromLog10(new double[]{log10PriorsByAC[0], MathUtils.log10sumLog10(log10PriorsByAC, 1)}, true);

        // TODO -- replace with more meaningful computation
        // TODO -- refactor this calculation into the ref calculation
        final Map<Allele, Double> log10pNonRefByAllele = new HashMap<Allele, Double>(allelesUsedInGenotyping.size());
        for ( int i = 0; i < subACOfMLE.length; i++ ) {
            final Allele allele = allelesUsedInGenotyping.get(i+1);
            final double log10PNonRef = getAlleleCountsOfMAP()[i] > 0 ? 0 : -10000; // TODO -- a total hack but in effect what the old behavior was
            log10pNonRefByAllele.put(allele, log10PNonRef);
        }

        return new AFCalcResult(subACOfMLE, nEvaluations, allelesUsedInGenotyping, log10Likelihoods, log10Priors, log10pNonRefByAllele);
    }

    // --------------------------------------------------------------------------------
    //
    // Protected mutational methods only for use within the calculation models themselves
    //
    // --------------------------------------------------------------------------------

    /**
     * Reset the data in this results object, so that it can be used in a subsequent AF calculation
     *
     * Resetting of the data is done by the calculation model itself, so shouldn't be done by callers any longer
     */
    protected void reset() {
        log10MLE = log10MAP = log10LikelihoodOfAFzero = log10PosteriorOfAFzero = VALUE_NOT_CALCULATED;
        for ( int i = 0; i < alleleCountsOfMLE.length; i++ ) {
            alleleCountsOfMLE[i] = 0;
            alleleCountsOfMAP[i] = 0;
        }
        currentLikelihoodsCacheIndex = 0;
        log10LikelihoodsMatrixSum = null;
        allelesUsedInGenotyping = null;
        nEvaluations = 0;
        Arrays.fill(log10LikelihoodsMatrixValues, Double.POSITIVE_INFINITY);
    }

    /**
     * Tell this result we used one more evaluation cycle
     */
    protected void incNEvaluations() {
        nEvaluations++;
    }

    protected void updateMLEifNeeded(final double log10LofK, final int[] alleleCountsForK) {
        addToLikelihoodsCache(log10LofK);

        if ( log10LofK > log10MLE ) {
            log10MLE = log10LofK;
            for ( int i = 0; i < alleleCountsForK.length; i++ )
                alleleCountsOfMLE[i] = alleleCountsForK[i];
        }
    }

    protected void updateMAPifNeeded(final double log10LofK, final int[] alleleCountsForK) {
        if ( log10LofK > log10MAP ) {
            log10MAP = log10LofK;
            for ( int i = 0; i < alleleCountsForK.length; i++ )
                alleleCountsOfMAP[i] = alleleCountsForK[i];
        }
    }

    private void addToLikelihoodsCache(final double log10LofK) {
        // add to the cache
        log10LikelihoodsMatrixValues[currentLikelihoodsCacheIndex++] = log10LofK;

        // if we've filled up the cache, then condense by summing up all of the values and placing the sum back into the first cell
        if ( currentLikelihoodsCacheIndex == LIKELIHOODS_CACHE_SIZE) {
            final double temporarySum = MathUtils.log10sumLog10(log10LikelihoodsMatrixValues, 0, currentLikelihoodsCacheIndex);
            Arrays.fill(log10LikelihoodsMatrixValues, Double.POSITIVE_INFINITY);
            log10LikelihoodsMatrixValues[0] = temporarySum;
            currentLikelihoodsCacheIndex = 1;
        }
    }

    protected void setLog10LikelihoodOfAFzero(final double log10LikelihoodOfAFzero) {
        this.log10LikelihoodOfAFzero = log10LikelihoodOfAFzero;
        if ( log10LikelihoodOfAFzero > log10MLE ) {
            log10MLE = log10LikelihoodOfAFzero;
            Arrays.fill(alleleCountsOfMLE, 0);
        }
    }

    protected void setLog10PosteriorOfAFzero(final double log10PosteriorOfAFzero) {
        this.log10PosteriorOfAFzero = log10PosteriorOfAFzero;
        if ( log10PosteriorOfAFzero > log10MAP ) {
            log10MAP = log10PosteriorOfAFzero;
            Arrays.fill(alleleCountsOfMAP, 0);
        }
    }

    protected void setAllelesUsedInGenotyping(List<Allele> allelesUsedInGenotyping) {
        if ( allelesUsedInGenotyping == null || allelesUsedInGenotyping.isEmpty() )
            throw new IllegalArgumentException("allelesUsedInGenotyping cannot be null or empty");

        this.allelesUsedInGenotyping = allelesUsedInGenotyping;
    }

    protected void setAClimits(int[] AClimits) {
        this.AClimits = AClimits;
    }
}