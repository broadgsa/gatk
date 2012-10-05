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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;


/**
 * Uses the Exact calculation of Heng Li
 */
abstract class ExactAFCalculation extends AlleleFrequencyCalculation {
    private final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6

    protected ExactAFCalculation(final UnifiedArgumentCollection UAC, final int nSamples, final Logger logger, final PrintStream verboseWriter) {
        super(UAC, nSamples, logger, verboseWriter);
    }

    protected ExactAFCalculation(final int nSamples, int maxAltAlleles, int maxAltAllelesForIndels, File exactCallsLog, Logger logger, PrintStream verboseWriter) {
        super(nSamples, maxAltAlleles, maxAltAllelesForIndels, exactCallsLog, logger, verboseWriter);
    }

    /**
     * Wrapper class that compares two likelihoods associated with two alleles
     */
    protected static final class LikelihoodSum implements Comparable<LikelihoodSum> {
        public double sum = 0.0;
        public Allele allele;

        public LikelihoodSum(Allele allele) { this.allele = allele; }

        public int compareTo(LikelihoodSum other) {
            final double diff = sum - other.sum;
            return ( diff < 0.0 ) ? 1 : (diff > 0.0 ) ? -1 : 0;
        }
    }

    /**
     * Unpack GenotypesContext into arraylist of doubel values
     * @param GLs            Input genotype context
     * @return               ArrayList of doubles corresponding to GL vectors
     */
    protected static ArrayList<double[]> getGLs(GenotypesContext GLs) {
        ArrayList<double[]> genotypeLikelihoods = new ArrayList<double[]>(GLs.size());

        genotypeLikelihoods.add(new double[]{0.0,0.0,0.0}); // dummy
        for ( Genotype sample : GLs.iterateInSampleNameOrder() ) {
            if ( sample.hasLikelihoods() ) {
                double[] gls = sample.getLikelihoods().getAsVector();

                if ( MathUtils.sum(gls) < UnifiedGenotyperEngine.SUM_GL_THRESH_NOCALL )
                    genotypeLikelihoods.add(gls);
            }
        }

        return genotypeLikelihoods;
    }

    /**
     * Computes the maximum ACs we need to consider for each alt allele
     *
     * Walks over the genotypes in VC, and computes for each alt allele the maximum
     * AC we need to consider in that alt allele dimension.  Does the calculation
     * based on the PLs in each genotype g, choosing to update the max AC for the
     * alt alleles corresponding to that PL.  Only takes the first lowest PL,
     * if there are multiple genotype configurations with the same PL value.  It
     * takes values in the order of the alt alleles.
     *
     * @param vc the variant context we will compute max alt alleles for
     * @return a vector of max alt alleles, indexed by alt allele, so result[0] is the AC of the
     *          first alt allele.
     */
    @Ensures("result != null")
    protected int[] computeMaxACs(final VariantContext vc) {
        final int[] maxACs = new int[vc.getNAlleles()-1];

        for ( final Genotype g : vc.getGenotypes() )
            updateMaxACs(g, maxACs);

        return maxACs;
    }

    /**
     * Update the maximum achievable allele counts in maxAC according to the PLs in g
     *
     * Selects the maximum genotype configuration from the PLs in g, and updates
     * the maxAC for this configure.  For example, if the lowest PL is for 0/1, updates
     * the maxAC for the alt allele 1 by 1.  If it's 1/1, update is 2.  Works for
     * many number of alt alleles (determined by length of maxACs).
     *
     * If the max PL occurs at 0/0, updates nothing
     * Note that this function greedily takes the first min PL, so that if 0/1 and 1/1 have
     * the same PL value, then updates the first one.
     *
     * Also, only will update 1 alt allele, so if 0/1 and 0/2 both have the same PL,
     * then only first one (1) will be updated
     *
     * @param g the genotype to update
     * @param maxACs the max allele count vector for alt alleles (starting at 0 => first alt allele)
     */
    @Requires({
            "g != null",
            "maxACs != null",
            "MathUtils.sum(maxACs) >= 0"})
    private void updateMaxACs(final Genotype g, final int[] maxACs) {
        final int[] PLs = g.getLikelihoods().getAsPLs();

        int minPLi = 0;
        int minPL = PLs[0];

        for ( int i = 0; i < PLs.length; i++ ) {
            if ( PLs[i] < minPL ) {
                minPL = PLs[i];
                minPLi = i;
            }
        }

        final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair pair = GenotypeLikelihoods.getAllelePair(minPLi);
        updateMaxACs(maxACs, pair.alleleIndex1);
        updateMaxACs(maxACs, pair.alleleIndex2);
    }

    /**
     * Simple helper.  Update max alt alleles maxACs according to the allele index (where 0 == ref)
     *
     * If alleleI == 0 => doesn't update anything
     * else maxACs[alleleI - 1]++
     *
     * @param maxACs array of max alt allele ACs
     * @param alleleI the index (relative to 0) to update a count of 1 in max alt alleles.
     */
    @Requires({
            "alleleI >= 0",
            "(alleleI - 1) < maxACs.length",
            "MathUtils.sum(maxACs) >= 0"})
    private void updateMaxACs(final int[] maxACs, final int alleleI) {
        if ( alleleI > 0 )
            maxACs[alleleI-1]++;
    }

    // -------------------------------------------------------------------------------------
    //
    // protected classes used to store exact model matrix columns
    //
    // -------------------------------------------------------------------------------------

    protected static final int HOM_REF_INDEX = 0;  // AA likelihoods are always first

    // a wrapper around the int array so that we can make it hashable
    protected static final class ExactACcounts {

        protected final int[] counts;
        private int hashcode = -1;

        public ExactACcounts(final int[] counts) {
            this.counts = counts;
        }

        public int[] getCounts() {
            return counts;
        }

        @Override
        public boolean equals(Object obj) {
            return (obj instanceof ExactACcounts) && Arrays.equals(counts, ((ExactACcounts) obj).counts);
        }

        @Override
        public int hashCode() {
            if ( hashcode == -1 )
                hashcode = Arrays.hashCode(counts);
            return hashcode;
        }

        @Override
        public String toString() {
            StringBuffer sb = new StringBuffer();
            sb.append(counts[0]);
            for ( int i = 1; i < counts.length; i++ ) {
                sb.append("/");
                sb.append(counts[i]);
            }
            return sb.toString();
        }
    }

    // This class represents a column in the Exact AC calculation matrix
    protected static final class ExactACset {

        // the counts of the various alternate alleles which this column represents
        final ExactACcounts ACcounts;

        // the column of the matrix
        final double[] log10Likelihoods;

        int sum = -1;

        public ExactACset(final int size, final ExactACcounts ACcounts) {
            this.ACcounts = ACcounts;
            log10Likelihoods = new double[size];
            Arrays.fill(log10Likelihoods, Double.NEGATIVE_INFINITY);
        }

        // sum of all the non-reference alleles
        public int getACsum() {
            if ( sum == -1 ) {
                sum = 0;
                for ( int count : ACcounts.getCounts() )
                    sum += count;
            }
            return sum;
        }

        public boolean equals(Object obj) {
            return (obj instanceof ExactACset) && ACcounts.equals(((ExactACset)obj).ACcounts);
        }
    }

    @Deprecated
    protected static final class OldMaxLikelihoodSeen {
        double maxLog10L = Double.NEGATIVE_INFINITY;
        ExactACcounts ACs = null;

        public OldMaxLikelihoodSeen() {}

        public void update(final double maxLog10L, final ExactACcounts ACs) {
            this.maxLog10L = maxLog10L;
            this.ACs = ACs;
        }

        // returns true iff all ACs in this object are less than or equal to their corresponding ACs in the provided set
        public boolean isLowerAC(final ExactACcounts otherACs) {
            final int[] myACcounts = this.ACs.getCounts();
            final int[] otherACcounts = otherACs.getCounts();

            for ( int i = 0; i < myACcounts.length; i++ ) {
                if ( myACcounts[i] > otherACcounts[i] )
                    return false;
            }
            return true;
        }
    }

    protected static final class MaxLikelihoodSeen {
        double maxLog10L = Double.NEGATIVE_INFINITY;
        final int[] maxACsToConsider;

        public MaxLikelihoodSeen(final int[] maxACsToConsider) {
            this.maxACsToConsider = maxACsToConsider;
        }

        /**
         * Update the maximum log10L seen, if log10LofKs is higher
         *
         * @param log10LofKs the likelihood of our current configuration state
         */
        public void update(final double log10LofKs) {
            if ( log10LofKs > maxLog10L )
                this.maxLog10L = log10LofKs;
        }

        /**
         * Is the likelihood of configuration K too low to consider, related to the
         * maximum likelihood seen already?
         *
         * @param log10LofK the log10 likelihood of the configuration we're considering analyzing
         * @return true if the configuration cannot meaningfully contribute to our likelihood sum
         */
        public boolean tooLowLikelihood(final double log10LofK) {
            return log10LofK < maxLog10L - MAX_LOG10_ERROR_TO_STOP_EARLY;
        }

        /**
         * Are all ACs in otherACs less than or equal to their corresponding ACs in the maxACsToConsider?
         *
         * @param otherACs the set of otherACs that we want to know if we should consider analyzing
         * @return true if otherACs is a state worth considering, or false otherwise
         */
        public boolean withinMaxACs(final ExactACcounts otherACs) {
            final int[] otherACcounts = otherACs.getCounts();

            for ( int i = 0; i < maxACsToConsider.length; i++ ) {
                // consider one more than the max AC to collect a bit more likelihood mass
                if ( otherACcounts[i] > maxACsToConsider[i] + 1 )
                    return false;
            }

            return true;
        }
    }
}