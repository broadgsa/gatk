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

import org.broadinstitute.sting.utils.MathUtils;

import java.util.Arrays;

public class DiploidSNPGenotypePriors implements GenotypePriors {
    // --------------------------------------------------------------------------------------------------------------
    //
    // Constants and static information
    //
    // --------------------------------------------------------------------------------------------------------------
    public static final double HUMAN_HETEROZYGOSITY = 1e-3;
    public static final double CEU_HETEROZYGOSITY = 1e-3;
    public static final double YRI_HETEROZYGOSITY = 1.0 / 850;

    /**
     * Default value of the prob of seeing a reference error.  Used to calculation the
     * chance of seeing a true B/C het when the reference is A, which we assume is the product
     * of the ref error rate and the het. Default value is Q60
     */
    public static final double PROB_OF_REFERENCE_ERROR = 1e-6;  // the reference is

    private final static double[] flatPriors = new double[DiploidGenotype.values().length];

    // --------------------------------------------------------------------------------------------------------------
    //
    // Diploid priors
    //
    // --------------------------------------------------------------------------------------------------------------
    private double[] priors = null;

    // todo -- fix me when this issue is resolved
    public static final boolean requirePriorSumToOne = false;

    /**
     * Create a new DiploidGenotypePriors object with flat priors for each diploid genotype
     */
    public DiploidSNPGenotypePriors() {
        priors = flatPriors.clone();
    }

    /**
     * Create a new GenotypeLikelihoods object with priors for a diploid with heterozygosity and reference
     * base ref
     *
     * @param ref
     * @param heterozygosity
     * @param probOfTriStateGenotype The prob of seeing a true B/C het when the reference is A
     */
    public DiploidSNPGenotypePriors(byte ref, double heterozygosity, double probOfTriStateGenotype) {
        priors = getReferencePolarizedPriors(ref, heterozygosity, probOfTriStateGenotype);
    }

    /**
     * Create a new Genotypelike Likelhoods's object with priors (in log10 space) for each of the DiploteGenotypes
     *
     * @param log10Priors
     */
    public DiploidSNPGenotypePriors(double[] log10Priors) {
        priors = log10Priors.clone();
    }

    /**
     * Returns an array of priors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return log10 prior as a double array
     */
    public double[] getPriors() {
        return priors;
    }

    /**
     * Returns the prior associated with DiploidGenotype g
     * @param g
     * @return log10 prior as a double
     */
    public double getPrior(DiploidGenotype g) {
        return getPriors()[g.ordinal()];
    }

    public double getHeterozygosity() { return HUMAN_HETEROZYGOSITY; }

    public boolean validate(boolean throwException) {
        try {
            if ( requirePriorSumToOne && MathUtils.compareDoubles(MathUtils.sumLog10(priors), 1.0) != 0 ) {
                throw new IllegalStateException(String.format("Priors don't sum to 1: sum=%f %s", MathUtils.sumLog10(priors), Arrays.toString(priors)));
            }

            for ( DiploidGenotype g : DiploidGenotype.values() ) {
                int i = g.ordinal();
                if ( ! MathUtils.wellFormedDouble(priors[i]) || ! MathUtils.isNegativeOrZero(priors[i]) ) {
                    String bad = String.format("Prior %f is badly formed %b", priors[i], MathUtils.isNegativeOrZero(priors[i]));
                    throw new IllegalStateException(String.format("At %s: %s", g.toString(), bad));
                }
            }
        } catch ( IllegalStateException e ) {
            if ( throwException )
                throw new RuntimeException(e);
            else
                return false;
        }

        return true;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Static functionality
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Returns homozygous-reference, heterozygous, and homozygous-non-ref probabilities given a heterozygosity
     * value, as elements 0, 1, and 2 of a double[], respectively
     *
     * @param h the heterozygosity [probability of a base being heterozygous]
     */
    @Deprecated
    public static double[] heterozygosity2DiploidProbabilities(double h) {
        double[] pdbls = new double[3];

        pdbls[0] = heterozygosity2HomRefProbability(h);
        pdbls[1] = heterozygosity2HetProbability(h);
        pdbls[2] = heterozygosity2HomVarProbability(h);
        return pdbls;
    }

    /**
     *
     * @param h
     * @return
     */
    public static double heterozygosity2HomRefProbability(double h) {
        if (MathUtils.isNegative(h)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        double v = 1.0 - (3.0 * h / 2.0);
        if (MathUtils.isNegative(v)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        return v;
    }

    public static double heterozygosity2HetProbability(double h) {
        if (MathUtils.isNegative(h)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        return h;
    }

    public static double heterozygosity2HomVarProbability(double h) {
        if (MathUtils.isNegative(h)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        return h / 2.0;
    }


    /**
     * Takes reference base, and three priors for hom-ref, het, hom-var, and fills in the priors vector
     * appropriately.
     *
     * Suppose A is the reference base, and we are given the probability of being hom-ref, het, and hom-var,
     * and that pTriSateGenotype is the true probability of observing reference A and a true genotype of B/C
     * then this sets the priors to:
     *
     * AA = hom-ref
     * AC = AG = AT = (het - pTriStateGenotype) / 3
     * CC = GG = TT = hom-var / 3
     * CG = CT = GT = pTriStateGenotype / 3
     *
     * So that we get:
     *
     * hom-ref + 3 * (het - pTriStateGenotype) / 3 + 3 * hom-var / 3 + 3 * pTriStateGenotype
     * hom-ref + het - pTriStateGenotype + hom-var + pTriStateGenotype
     * hom-ref + het + hom-var
     * = 1
     *
     * @param ref
     * @param heterozyosity
     * @param pRefError
     */
    public static double[] getReferencePolarizedPriors(byte ref, double heterozyosity, double pRefError ) {
        if ( ! MathUtils.isBounded(pRefError, 0.0, 0.01) ) {
            throw new RuntimeException(String.format("BUG: p Reference error is out of bounds (0.0 - 0.01) is allow range %f", pRefError));
        }

        double pTriStateGenotype = heterozyosity * pRefError;
//        if ( pTriStateGenotype >= heterozyosity ) {
//            throw new RuntimeException(String.format("p Tristate genotype %f is greater than the heterozygosity %f", pTriStateGenotype, heterozyosity));
//        }

        double pHomRef = heterozygosity2HomRefProbability(heterozyosity);
        double pHet    = heterozygosity2HetProbability(heterozyosity);
        double pHomVar = heterozygosity2HomVarProbability(heterozyosity);

        if (MathUtils.compareDoubles(pHomRef + pHet + pHomVar, 1.0) != 0) {
            throw new RuntimeException(String.format("BUG: Prior probabilities don't sum to one => %f, %f, %f", pHomRef, pHet, pHomVar));
        }

        double[] priors = new double[DiploidGenotype.values().length];

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double POfG;

            final double nOnRefHets = 3;
            final double nOffRefHets = 3;
            final double nHomVars = 3;

            if ( g.isHomRef(ref) )      { POfG = pHomRef; }
            else if ( g.isHomVar(ref) ) { POfG = pHomVar / nHomVars; }
            else if ( g.isHetRef(ref) ) { POfG = (pHet - pTriStateGenotype ) / nOnRefHets; }
            else                        { POfG = pTriStateGenotype / nOffRefHets; }

            priors[g.ordinal()] = Math.log10(POfG);
        }

        return priors;
    }

    static {
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            flatPriors[g.ordinal()] = Math.log10(1.0 / DiploidGenotype.values().length);
        }
    }
}
