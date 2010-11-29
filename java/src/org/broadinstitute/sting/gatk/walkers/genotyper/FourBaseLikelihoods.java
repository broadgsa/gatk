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

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.*;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

/**
 * Stable, error checking version of strict 4 base likelihoods.  Useful for calculating the likelihoods, priors,
 * and posteriors given a pile of bases and quality scores (in conjuncion with GenotypeLikelihoods)
 *
 * Suppose we have bases b1, b2, ..., bN with qualities scores q1, q2, ..., qN.  This object
 * calculates:
 *
 * P(b | D) = P(b) * P(D | b)
 *
 * where
 *
 * P(D | b) = sum_i log10 P(bi | b)
 *
 * and
 *
 * P(bi | b) = 1 - P(error | q1) if bi = b
 *           = P(error | q1) / 3 if bi != b
 *
 *
 */
public abstract class FourBaseLikelihoods implements Cloneable {

    protected boolean enableCacheFlag = true;

    //
    // The fundamental data array associated with 4-base likelihoods
    //
    protected double[] log10Likelihoods = null;


    /**
     * If true, lots of output will be generated about the Likelihoods at each site
     */
    private boolean verbose = false;

    /**
     * Bases with Q scores below this threshold aren't included in the Likelihood calculation
     */
    private int minQScoreToInclude = 0;

    /**
     * Create a new FourBaseLikelihoods object
     */
    public FourBaseLikelihoods() {
        log10Likelihoods = zeros.clone();          // Likelihoods are all zeros
    }

    /**
     * Cloning of the object
     * @return clone
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        FourBaseLikelihoods c = (FourBaseLikelihoods)super.clone();
        c.log10Likelihoods = log10Likelihoods.clone();
        return c;
    }

    public void setVerbose(boolean v) {
        verbose = v;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public int getMinQScoreToInclude() {
        return minQScoreToInclude;
    }

    public void setMinQScoreToInclude(int minQScoreToInclude) {
        this.minQScoreToInclude = minQScoreToInclude;
    }

    /**
     * Returns an array of log10 likelihoods for each base, indexed by BaseUtils.BASES.ordinal values()
     * @return probs
     */
    public double[] getLog10Likelihoods() {
        return log10Likelihoods;
    }

    /**
     * Returns the likelihood associated with a base
     * @param base base
     * @return log10 likelihood as a double
     */
    public double getLog10Likelihood(byte base) {
        return getLog10Likelihood(BaseUtils.simpleBaseToBaseIndex(base));
    }

    public double getLog10Likelihood(int baseIndex) {
        return (baseIndex < 0 ? 0.0 : getLog10Likelihoods()[baseIndex]);
    }

    /**
     * Returns an array of likelihoods for each base, indexed by BaseUtils.BASES.ordinal values()
     * @return probs
     */
    public double[] getLikelihoods() {
        double[] probs = new double[4];
        for (int i = 0; i < 4; i++)
            probs[i] = Math.pow(10, log10Likelihoods[i]);
        return probs;
    }

    /**
     * Returns the likelihoods associated with a base
     * @param base base
     * @return likelihoods as a double
     */
    public double getLikelihood(byte base) {
        return getLikelihood(BaseUtils.simpleBaseToBaseIndex(base));
    }

    public double getLikelihood(int baseIndex) {
        return (baseIndex < 0 ? 0.0 : Math.pow(10, log10Likelihoods[baseIndex]));
    }


    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // add() -- the heart of
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Updates likelihoods and posteriors to reflect an additional observation of observedBase with
     * qualityScore.
     *
     * @param observedBase observed base
     * @param qualityScore base quality
     * @param read         SAM read
     * @param offset       offset on read
     * @return 1 if the base was considered good enough to add to the likelihoods (not Q0 or 'N', for example)
     */
    public int add(byte observedBase, byte qualityScore, SAMRecord read, int offset) {
        FourBaseLikelihoods fbl = computeLog10Likelihoods(observedBase, qualityScore, read, offset);
        if ( fbl == null )
            return 0;

        for ( byte base : BaseUtils.BASES ) {
            double likelihood = fbl.getLikelihood(base);
            log10Likelihoods[BaseUtils.simpleBaseToBaseIndex(base)] += likelihood;
        }

        if ( verbose ) {
            for ( byte base : BaseUtils.BASES ) { System.out.printf("%s\t", (char)base); }
            System.out.println();
            for ( byte base : BaseUtils.BASES ) { System.out.printf("%.2f\t", log10Likelihoods[BaseUtils.simpleBaseToBaseIndex(base)]); }
            System.out.println();
        }

        return 1;
    }

    /**
     * Updates likelihoods and posteriors to reflect an additional observation of observedBase with
     * qualityScore.
     *
     * @param observedBase observed base
     * @param qualityScore base quality
     * @param read         SAM read
     * @param offset       offset on read
     * @return likelihoods for this observation or null if the base was not considered good enough to add to the likelihoods (Q0 or 'N', for example)
     */
    public FourBaseLikelihoods computeLog10Likelihoods(byte observedBase, byte qualityScore, SAMRecord read, int offset) {
        if ( badBase(observedBase) ) {
            return null;
        }

        try {
            if ( qualityScore > getMinQScoreToInclude() ) {

                FourBaseLikelihoods fbl = (FourBaseLikelihoods)this.clone();
                fbl.log10Likelihoods = zeros.clone();

                for ( byte base : BaseUtils.BASES ) {
                    double likelihood = log10PofObservingBaseGivenChromosome(observedBase, base, qualityScore, read, offset);

                    if ( verbose ) {
                        boolean fwdStrand = ! read.getReadNegativeStrandFlag();
                        System.out.printf("  L(%c | b=%s, Q=%d, S=%s) = %f / %f%n",
                                observedBase, base, qualityScore, fwdStrand ? "+" : "-", pow(10,likelihood) * 100, likelihood);
                    }

                    fbl.log10Likelihoods[BaseUtils.simpleBaseToBaseIndex(base)] = likelihood;
                }
                
                return fbl;
            }
        } catch ( CloneNotSupportedException e ) {
            throw new RuntimeException(e);
        }

        return null;
    }


    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // helper routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Returns true when the observedBase is considered bad and shouldn't be processed by this object.  A base
     * is considered bad if:
     *
     *   Criterion 1: observed base isn't a A,C,T,G or lower case equivalent
     *
     * @param observedBase observed base
     * @return true if the base is a bad base
     */
    private boolean badBase(byte observedBase) {
        return BaseUtils.simpleBaseToBaseIndex(observedBase) == -1;
    }

    /**
     * Return a string representation of this object in a moderately usable form
     *
     * @return string representation
     */
    public String toString() {
        double sum = 0;
        StringBuilder s = new StringBuilder();
        for ( byte base : BaseUtils.BASES ) {
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);
            s.append(String.format("%s %.10f ", base, log10Likelihoods[baseIndex]));
			sum += Math.pow(10, log10Likelihoods[baseIndex]);
        }
		s.append(String.format(" %f", sum));
        return s.toString();
    }

    // in general, we don't care about the platform index; EmpiricalSubstitutionProbabilities overloads this
    public int getReadSequencerPlatformIndex( SAMRecord read ) {
        return 0;
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // Validation routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    public boolean validate() {
        return validate(true);
    }

    public boolean validate(boolean throwException) {
        try {

            for ( byte base : BaseUtils.BASES ) {

                int i = BaseUtils.simpleBaseToBaseIndex(base);
                if ( ! MathUtils.wellFormedDouble(log10Likelihoods[i]) || ! MathUtils.isNegativeOrZero(log10Likelihoods[i]) ) {
                    String bad = String.format("Likelihood %f is badly formed", log10Likelihoods[i]);
                    throw new IllegalStateException(String.format("At %s: %s", base, bad));
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

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // Hearty math calculations follow
    //
    //    -- these should not be messed with unless you know what you are doing
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     *
     * @param observedBase observed base
     * @param chromBase    target base
     * @param qual         base quality
     * @param read         SAM read
     * @param offset       offset on read
     * @return log10 likelihood
     */

    protected double log10PofObservingBaseGivenChromosome(byte observedBase, byte chromBase, byte qual, SAMRecord read, int offset) {
        if (qual == 0) { // zero quals are wrong
            throw new RuntimeException(String.format("Unexpected Q0 base discovered in log10PofObservingBaseGivenChromosome: %c %s %d at %d in %s",
                    observedBase, chromBase, qual, offset, read));
        }

        double logP;

        if ( observedBase == chromBase ) {
            // the base is consistent with the chromosome -- it's 1 - e
            //logP = oneMinusData[qual];
            double e = pow(10, (qual / -10.0));
            logP = log10(1.0 - e);
        } else {
            // the base is inconsistent with the chromosome -- it's e * P(chromBase | observedBase is an error)
            logP = qual / -10.0 + log10PofTrueBaseGivenMiscall(observedBase, chromBase, read, offset);
        }

        //System.out.printf("%c %c %d => %f%n", observedBase, chromBase, qual, logP);
        return logP;
    }

    /**
     * Must be overridden by concrete subclasses
     *
     * @param observedBase observed base
     * @param chromBase    target base
     * @param read         SAM read
     * @param offset       offset on read
    * @return log10 likelihood
     */
    protected abstract double log10PofTrueBaseGivenMiscall(byte observedBase, byte chromBase, SAMRecord read, int offset);

    //
    // Constant static data
    //
    private final static double[] zeros = new double[BaseUtils.BASES.length];

    static {
        for ( byte base : BaseUtils.BASES ) {
            zeros[BaseUtils.simpleBaseToBaseIndex(base)] = 0.0;
        }
    }
}