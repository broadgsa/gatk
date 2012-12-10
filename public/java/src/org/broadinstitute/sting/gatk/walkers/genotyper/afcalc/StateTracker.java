package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Keeps track of the state information during the exact model AF calculation.
 *
 * Tracks things like the MLE and MAP AC values, their corresponding likelhood and posterior
 * values, the likelihood of the AF == 0 state, and the number of evaluations needed
 * by the calculation to compute the P(AF == 0)
 */
final class StateTracker {
    protected static final double VALUE_NOT_CALCULATED = Double.NEGATIVE_INFINITY;
    protected final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6

    /**
     * These variables are intended to contain the MLE and MAP (and their corresponding allele counts)
     * of the site over all alternate alleles
     */
    protected double log10MLE;
    protected double log10MAP;

    /**
     * Returns a vector with maxAltAlleles values containing AC values at the MLE
     *
     * The values of the ACs for this call are stored in the getAllelesUsedInGenotyping order,
     * starting from index 0 (i.e., the first alt allele is at 0).  The vector is always
     * maxAltAlleles in length, and so only the first getAllelesUsedInGenotyping.size() - 1 values
     * are meaningful.
     */
    private final int[] alleleCountsOfMLE;
    private final int[] alleleCountsOfMAP;

    /**
     * A vector of log10 likelihood values seen, for future summation.  When the size of the
     * vector is exceeed -- because we've pushed more posteriors than there's space to hold
     * -- we simply sum up the existing values, make that the first value, and continue.
     */
    private final double[] log10LikelihoodsForAFGt0 = new double[LIKELIHOODS_CACHE_SIZE];
    private static final int LIKELIHOODS_CACHE_SIZE = 5000;
    private int log10LikelihoodsForAFGt0CacheIndex = 0;

    /**
     * The actual sum of the likelihoods.  Null if the sum hasn't been computed yet
     */
    protected Double log10LikelihoodsForAFGt0Sum = null;

    /**
     * Contains the likelihood for the site's being monomorphic (i.e. AF=0 for all alternate alleles)
     */
    private double log10LikelihoodOfAFzero = 0.0;

    /**
     * The number of evaluates we've gone through in the AFCalc
     */
    private int nEvaluations = 0;

    /**
     * The list of alleles actually used in computing the AF
     */
    private List<Allele> allelesUsedInGenotyping = null;

    /**
     * Create a results object capability of storing results for calls with up to maxAltAlleles
     *
     * @param maxAltAlleles an integer >= 1
     */
    public StateTracker(final int maxAltAlleles) {
        if ( maxAltAlleles < 1 ) throw new IllegalArgumentException("maxAltAlleles must be >= 0, saw " + maxAltAlleles);

        alleleCountsOfMLE = new int[maxAltAlleles];
        alleleCountsOfMAP = new int[maxAltAlleles];

        reset();
    }

    /**
     * Is the likelihood of configuration K too low to consider, related to the
     * maximum likelihood seen already?
     *
     * @param log10LofK the log10 likelihood of the configuration we're considering analyzing
     * @return true if the configuration cannot meaningfully contribute to our likelihood sum
     */
    private boolean tooLowLikelihood(final double log10LofK) {
        return log10LofK < log10MLE - MAX_LOG10_ERROR_TO_STOP_EARLY;
    }

    /**
     * @return true iff all ACs in this object are less than or equal to their corresponding ACs in the provided set
     */
    private boolean isLowerAC(final ExactACcounts otherACs) {
        final int[] otherACcounts = otherACs.getCounts();

        for ( int i = 0; i < otherACcounts.length; i++ ) {
            if ( alleleCountsOfMLE[i] > otherACcounts[i] )
                return false;
        }

        return true;
    }

    /**
     * Should we stop exploring paths from ACs, given it's log10LofK
     *
     * @param log10LofK the log10LofK of these ACs
     * @param ACs the ACs of this state
     * @return return true if there's no reason to continue with subpaths of AC, or false otherwise
     */
    protected boolean abort( final double log10LofK, final ExactACcounts ACs, final boolean enforceLowerACs ) {
        return tooLowLikelihood(log10LofK) && (!enforceLowerACs || isLowerAC(ACs));
    }

    @Ensures("result != null")
    protected int[] getAlleleCountsOfMAP() {
        return alleleCountsOfMAP;
    }

    @Ensures("result >= 0")
    protected int getnEvaluations() {
        return nEvaluations;
    }

    /**
     * @return the likelihoods summed across all AC values for AC > 0
     */
    private double getLog10LikelihoodOfAFNotZero() {
        if ( log10LikelihoodsForAFGt0Sum == null ) {
            if ( log10LikelihoodsForAFGt0CacheIndex == 0 ) // there's nothing to sum up, so make the sum equal to the smallest thing we have
                log10LikelihoodsForAFGt0Sum = MathUtils.LOG10_P_OF_ZERO;
            else
                log10LikelihoodsForAFGt0Sum = MathUtils.log10sumLog10(log10LikelihoodsForAFGt0, 0, log10LikelihoodsForAFGt0CacheIndex);
        }
        return log10LikelihoodsForAFGt0Sum;
    }

    /**
     * @return the log10 likelihood of AF == 0
     */
    private double getLog10LikelihoodOfAFzero() {
        return log10LikelihoodOfAFzero;
    }

    /**
     * Convert this state to an corresponding AFCalcResult.
     *
     * Assumes that the values in this state have been filled in with meaningful values during the calculation.
     * For example, that the allelesUsedInGenotyping has been set, that the alleleCountsOfMLE contains meaningful
     * values, etc.
     *
     * @param log10PriorsByAC the priors by AC
     *
     * @return an AFCalcResult summarizing the final results of this calculation
     */
    @Requires("allelesUsedInGenotyping != null")
    protected AFCalcResult toAFCalcResult(final double[] log10PriorsByAC) {
        final int [] subACOfMLE = Arrays.copyOf(alleleCountsOfMLE, allelesUsedInGenotyping.size() - 1);
        final double[] log10Likelihoods = MathUtils.normalizeFromLog10(new double[]{getLog10LikelihoodOfAFzero(), getLog10LikelihoodOfAFNotZero()}, true);
        final double[] log10Priors = MathUtils.normalizeFromLog10(new double[]{log10PriorsByAC[0], MathUtils.log10sumLog10(log10PriorsByAC, 1)}, true);

        final Map<Allele, Double> log10pRefByAllele = new HashMap<Allele, Double>(allelesUsedInGenotyping.size());
        for ( int i = 0; i < subACOfMLE.length; i++ ) {
            final Allele allele = allelesUsedInGenotyping.get(i+1);
            final double log10PRef = alleleCountsOfMAP[i] > 0 ? -10000 : 0; // TODO -- a total hack but in effect what the old behavior was
            log10pRefByAllele.put(allele, log10PRef);
        }

        return new AFCalcResult(subACOfMLE, nEvaluations, allelesUsedInGenotyping, log10Likelihoods, log10Priors, log10pRefByAllele);
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
        log10MLE = log10MAP = log10LikelihoodOfAFzero = VALUE_NOT_CALCULATED;
        log10LikelihoodsForAFGt0CacheIndex = 0;
        log10LikelihoodsForAFGt0Sum = null;
        allelesUsedInGenotyping = null;
        nEvaluations = 0;
        Arrays.fill(alleleCountsOfMLE, 0);
        Arrays.fill(alleleCountsOfMAP, 0);
        Arrays.fill(log10LikelihoodsForAFGt0, Double.POSITIVE_INFINITY);
    }

    /**
     * Tell this result we used one more evaluation cycle
     */
    protected void incNEvaluations() {
        nEvaluations++;
    }

    /**
     * Update the maximum log10 likelihoods seen, if log10LofKs is higher, and the corresponding ACs of this state
     *
     * @param log10LofK the likelihood of our current configuration state, cannot be the 0 state
     * @param alleleCountsForK the allele counts for this state
     */
    @Requires({"alleleCountsForK != null", "MathUtils.sum(alleleCountsForK) >= 0"})
    @Ensures("log10MLE == Math.max(log10LofK, log10MLE)")
    protected void updateMLEifNeeded(final double log10LofK, final int[] alleleCountsForK) {
        addToLikelihoodsCache(log10LofK);

        if ( log10LofK > log10MLE ) {
            log10MLE = log10LofK;
            System.arraycopy(alleleCountsForK, 0, alleleCountsOfMLE, 0, alleleCountsForK.length);
        }
    }

    /**
     * Update the maximum log10 posterior seen, if log10PofKs is higher, and the corresponding ACs of this state
     *
     * @param log10PofK the posterior of our current configuration state
     * @param alleleCountsForK the allele counts for this state
     */
    @Requires({"alleleCountsForK != null", "MathUtils.sum(alleleCountsForK) >= 0"})
    @Ensures("log10MAP == Math.max(log10PofK, log10MAP)")
    protected void updateMAPifNeeded(final double log10PofK, final int[] alleleCountsForK) {
        if ( log10PofK > log10MAP ) {
            log10MAP = log10PofK;
            System.arraycopy(alleleCountsForK, 0, alleleCountsOfMAP, 0, alleleCountsForK.length);
        }
    }

    private void addToLikelihoodsCache(final double log10LofK) {
        // add to the cache
        log10LikelihoodsForAFGt0[log10LikelihoodsForAFGt0CacheIndex++] = log10LofK;

        // if we've filled up the cache, then condense by summing up all of the values and placing the sum back into the first cell
        if ( log10LikelihoodsForAFGt0CacheIndex == LIKELIHOODS_CACHE_SIZE) {
            final double temporarySum = MathUtils.log10sumLog10(log10LikelihoodsForAFGt0, 0, log10LikelihoodsForAFGt0CacheIndex);
            Arrays.fill(log10LikelihoodsForAFGt0, Double.POSITIVE_INFINITY);
            log10LikelihoodsForAFGt0[0] = temporarySum;
            log10LikelihoodsForAFGt0CacheIndex = 1;
        }
    }

    protected void setLog10LikelihoodOfAFzero(final double log10LikelihoodOfAFzero) {
        this.log10LikelihoodOfAFzero = log10LikelihoodOfAFzero;
        if ( log10LikelihoodOfAFzero > log10MLE ) {
            log10MLE = log10LikelihoodOfAFzero;
            Arrays.fill(alleleCountsOfMLE, 0);
        }
    }

    @Requires({"MathUtils.goodLog10Probability(log10PosteriorOfAFzero)"})
    protected void setLog10PosteriorOfAFzero(final double log10PosteriorOfAFzero) {
        if ( log10PosteriorOfAFzero > log10MAP ) {
            log10MAP = log10PosteriorOfAFzero;
            Arrays.fill(alleleCountsOfMAP, 0);
        }
    }

    /**
     * Set the list of alleles used in genotyping
     *
     * @param allelesUsedInGenotyping the list of alleles, where the first allele is reference
     */
    @Requires({"allelesUsedInGenotyping != null", "allelesUsedInGenotyping.size() > 1"})
    protected void setAllelesUsedInGenotyping(List<Allele> allelesUsedInGenotyping) {
        if ( allelesUsedInGenotyping == null || allelesUsedInGenotyping.isEmpty() )
            throw new IllegalArgumentException("allelesUsedInGenotyping cannot be null or empty");
        if ( allelesUsedInGenotyping.get(0).isNonReference() )
            throw new IllegalArgumentException("The first element of allelesUsedInGenotyping must be the reference allele");

        this.allelesUsedInGenotyping = allelesUsedInGenotyping;
    }
}
