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

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;


/**
 * Uses the Exact calculation of Heng Li
 */
abstract class ExactAFCalculation extends AlleleFrequencyCalculation {
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

    protected static final class MaxLikelihoodSeen {
        double maxLog10L = Double.NEGATIVE_INFINITY;
        ExactACcounts ACs = null;

        public MaxLikelihoodSeen() {}

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
}