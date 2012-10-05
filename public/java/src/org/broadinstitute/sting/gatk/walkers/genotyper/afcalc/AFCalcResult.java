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
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Dec 14, 2011
 *
 * Useful helper class to communicate the results of the allele frequency calculation
 *
 * TODO -- WHAT IS THE CONTRACT ON MAP AC AND P NON REF?
 */
public class AFCalcResult {
    // These variables are intended to contain the MLE and MAP (and their corresponding allele counts) of the site over all alternate alleles
    private double log10MLE;
    private double log10MAP;
    private final int[] alleleCountsOfMLE;
    private final int[] alleleCountsOfMAP;

    // The posteriors seen, not including that of AF=0
    private static final int POSTERIORS_CACHE_SIZE = 5000;
    private final double[] log10PosteriorMatrixValues = new double[POSTERIORS_CACHE_SIZE];
    private int currentPosteriorsCacheIndex = 0;
    private Double log10PosteriorMatrixSum = null;

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
    public AFCalcResult(final int maxAltAlleles) {
        if ( maxAltAlleles < 1 ) throw new IllegalArgumentException("maxAltAlleles must be >= 0, saw " + maxAltAlleles);

        alleleCountsOfMLE = new int[maxAltAlleles];
        alleleCountsOfMAP = new int[maxAltAlleles];

        reset();
    }

    /**
     * Get the log10 value of the probability mass at the MLE
     *
     * @return a log10 prob
     */
    @Ensures("goodLog10Value(result)")
    public double getLog10MLE() {
        return log10MLE;
    }

    /**
     * Get the log10 value of the probability mass at the max. a posterior (MAP)
     *
     * @return a log10 prob
     */
    @Ensures("goodLog10Value(result)")
    public double getLog10MAP() {
        return log10MAP;
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
     * Returns the number of cycles used to evaluate the pNonRef for this AF calculation
     *
     * @return the number of evaluations required to produce the answer for this AF calculation
     */
    public int getnEvaluations() {
        return nEvaluations;
    }

    /**
     * TODO -- eric what is this supposed to return?  my unit tests don't do what I think they should
     *
     * @return
     */
    public double getLog10PosteriorsMatrixSumWithoutAFzero() {
        if ( log10PosteriorMatrixSum == null ) {
            log10PosteriorMatrixSum = MathUtils.log10sumLog10(log10PosteriorMatrixValues, 0, currentPosteriorsCacheIndex);
        }
        return log10PosteriorMatrixSum;
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

    /**
     * Get the list of alleles actually used in genotyping.
     *
     * Due to computational / implementation constraints this may be smaller than
     * the actual list of alleles requested
     *
     * @return a non-empty list of alleles used during genotyping
     */
    @Ensures({"result != null", "! result.isEmpty()"})
    public List<Allele> getAllelesUsedInGenotyping() {
        if ( allelesUsedInGenotyping == null )
            throw new IllegalStateException("allelesUsedInGenotyping requested but not yet set");

        return allelesUsedInGenotyping;
    }

    /**
     * Get the normalized -- across all AFs -- of AC == 0, NOT LOG10
     * @return
     */
    // TODO -- this ensure cannot be enabled right now because the log10 inputs can be infinity, etc.
    // TODO -- we should own these values in a more meaningful way and return good values in the case
    // TODO -- where this happens, or instead thrown an error and have a function to say "was this calculation successful
//    @Ensures({"result >= 0.0", "result <= 1.0"})
    public double getNormalizedPosteriorOfAFzero() {
        return getNormalizedPosteriors()[0];
    }

    /**
     * Get the normalized -- across all AFs -- of AC > 0, NOT LOG10
     * @return
     */
    // TODO -- this ensure cannot be enabled right now because the log10 inputs can be infinity, etc.
    // TODO -- we should own these values in a more meaningful way and return good values in the case
    // TODO -- where this happens, or instead thrown an error and have a function to say "was this calculation successful
    //@Ensures({"result >= 0.0", "result <= 1.0"})
    public double getNormalizedPosteriorOfAFGTZero() {
        return getNormalizedPosteriors()[1];
    }

    private double[] getNormalizedPosteriors() {
        final double[] posteriors = new double[]{ getLog10PosteriorOfAFzero(), getLog10PosteriorsMatrixSumWithoutAFzero() };
        return MathUtils.normalizeFromLog10(posteriors);
    }

    public int[] getAClimits() {
        return AClimits;
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
        log10MLE = log10MAP = log10LikelihoodOfAFzero = log10PosteriorOfAFzero = AFCalc.VALUE_NOT_CALCULATED;
        for ( int i = 0; i < alleleCountsOfMLE.length; i++ ) {
            alleleCountsOfMLE[i] = 0;
            alleleCountsOfMAP[i] = 0;
        }
        currentPosteriorsCacheIndex = 0;
        log10PosteriorMatrixSum = null;
        allelesUsedInGenotyping = null;
    }

    /**
     * Tell this result we used one more evaluation cycle
     */
    protected void incNEvaluations() {
        nEvaluations++;
    }

    protected void updateMLEifNeeded(final double log10LofK, final int[] alleleCountsForK) {
        if ( log10LofK > log10MLE ) {
            log10MLE = log10LofK;
            for ( int i = 0; i < alleleCountsForK.length; i++ )
                alleleCountsOfMLE[i] = alleleCountsForK[i];
        }
    }

    protected void updateMAPifNeeded(final double log10LofK, final int[] alleleCountsForK) {
        addToPosteriorsCache(log10LofK);

        if ( log10LofK > log10MAP ) {
            log10MAP = log10LofK;
            for ( int i = 0; i < alleleCountsForK.length; i++ )
                alleleCountsOfMAP[i] = alleleCountsForK[i];
        }
    }

    private void addToPosteriorsCache(final double log10LofK) {
        // add to the cache
        log10PosteriorMatrixValues[currentPosteriorsCacheIndex++] = log10LofK;

        // if we've filled up the cache, then condense by summing up all of the values and placing the sum back into the first cell
        if ( currentPosteriorsCacheIndex == POSTERIORS_CACHE_SIZE ) {
            final double temporarySum = MathUtils.log10sumLog10(log10PosteriorMatrixValues, 0, currentPosteriorsCacheIndex);
            log10PosteriorMatrixValues[0] = temporarySum;
            currentPosteriorsCacheIndex = 1;
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

    private static boolean goodLog10Value(final double result) {
        return result <= 0.0 || Double.isInfinite(result) || Double.isNaN(result);
    }

    protected void setAClimits(int[] AClimits) {
        this.AClimits = AClimits;
    }
}