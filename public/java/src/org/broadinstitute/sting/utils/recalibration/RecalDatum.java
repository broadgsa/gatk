package org.broadinstitute.sting.utils.recalibration;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

import com.google.java.contract.Ensures;
import com.google.java.contract.Invariant;
import com.google.java.contract.Requires;
import org.apache.commons.math.optimization.fitting.GaussianFunction;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;


/**
 * An individual piece of recalibration data. Each bin counts up the number of observations and the number
 * of reference mismatches seen for that combination of covariates.
 *
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 3, 2009
 */
@Invariant({
        "estimatedQReported >= 0.0",
        "! Double.isNaN(estimatedQReported)",
        "! Double.isInfinite(estimatedQReported)",
        "empiricalQuality >= 0.0 || empiricalQuality == UNINITIALIZED",
        "! Double.isNaN(empiricalQuality)",
        "! Double.isInfinite(empiricalQuality)",
        "numObservations >= 0",
        "numMismatches >= 0",
        "numMismatches <= numObservations"
})
public class RecalDatum {
    private static final double UNINITIALIZED = -1.0;

    /**
     * estimated reported quality score based on combined data's individual q-reporteds and number of observations
     */
    private double estimatedQReported;

    /**
     * the empirical quality for datums that have been collapsed together (by read group and reported quality, for example)
     */
    private double empiricalQuality;

    /**
     * number of bases seen in total
     */
    private double numObservations;

    /**
     * number of bases seen that didn't match the reference
     */
    private double numMismatches;

    /**
     * used when calculating empirical qualities to avoid division by zero
     */
    private static final int SMOOTHING_CONSTANT = 1;

    //---------------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Create a new RecalDatum with given observation and mismatch counts, and an reported quality
     *
     * @param _numObservations
     * @param _numMismatches
     * @param reportedQuality
     */
    public RecalDatum(final double _numObservations, final double _numMismatches, final byte reportedQuality) {
        if ( _numObservations < 0 ) throw new IllegalArgumentException("numObservations < 0");
        if ( _numMismatches < 0 ) throw new IllegalArgumentException("numMismatches < 0");
        if ( reportedQuality < 0 ) throw new IllegalArgumentException("reportedQuality < 0");

        numObservations = _numObservations;
        numMismatches = _numMismatches;
        estimatedQReported = reportedQuality;
        empiricalQuality = UNINITIALIZED;
    }

    /**
     * Copy copy into this recal datum, overwriting all of this objects data
     * @param copy
     */
    public RecalDatum(final RecalDatum copy) {
        this.numObservations = copy.getNumObservations();
        this.numMismatches = copy.getNumMismatches();
        this.estimatedQReported = copy.estimatedQReported;
        this.empiricalQuality = copy.empiricalQuality;
    }

    /**
     * Add in all of the data from other into this object, updating the reported quality from the expected
     * error rate implied by the two reported qualities
     *
     * @param other
     */
    public synchronized void combine(final RecalDatum other) {
        final double sumErrors = this.calcExpectedErrors() + other.calcExpectedErrors();
        increment(other.getNumObservations(), other.getNumMismatches());
        estimatedQReported = -10 * Math.log10(sumErrors / getNumObservations());
        empiricalQuality = UNINITIALIZED;
    }

    public synchronized void setEstimatedQReported(final double estimatedQReported) {
        if ( estimatedQReported < 0 ) throw new IllegalArgumentException("estimatedQReported < 0");
        if ( Double.isInfinite(estimatedQReported) ) throw new IllegalArgumentException("estimatedQReported is infinite");
        if ( Double.isNaN(estimatedQReported) ) throw new IllegalArgumentException("estimatedQReported is NaN");

        this.estimatedQReported = estimatedQReported;
    }

    public final double getEstimatedQReported() {
        return estimatedQReported;
    }
    public final byte getEstimatedQReportedAsByte() {
        return (byte)(int)(Math.round(getEstimatedQReported()));
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // Empirical quality score -- derived from the num mismatches and observations
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Returns the error rate (in real space) of this interval, or 0 if there are no observations
     * @return the empirical error rate ~= N errors / N obs
     */
    @Ensures("result >= 0.0")
    public double getEmpiricalErrorRate() {
        if ( numObservations == 0 )
            return 0.0;
        else {
            // cache the value so we don't call log over and over again
            final double doubleMismatches = numMismatches + SMOOTHING_CONSTANT;
            // smoothing is one error and one non-error observation, for example
            final double doubleObservations = numObservations + SMOOTHING_CONSTANT + SMOOTHING_CONSTANT;
            return doubleMismatches / doubleObservations;
        }
    }

    public synchronized void setEmpiricalQuality(final double empiricalQuality) {
        if ( empiricalQuality < 0 ) throw new IllegalArgumentException("empiricalQuality < 0");
        if ( Double.isInfinite(empiricalQuality) ) throw new IllegalArgumentException("empiricalQuality is infinite");
        if ( Double.isNaN(empiricalQuality) ) throw new IllegalArgumentException("empiricalQuality is NaN");

        this.empiricalQuality = empiricalQuality;
    }

    public final double getEmpiricalQuality() {
        if (empiricalQuality == UNINITIALIZED)
            calcEmpiricalQuality();
        return empiricalQuality;
    }

    public final byte getEmpiricalQualityAsByte() {
        return (byte)(Math.round(getEmpiricalQuality()));
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // toString methods
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public String toString() {
        return String.format("%.2f,%.2f,%.2f", getNumObservations(), getNumMismatches(), getEmpiricalQuality());
    }

    public String stringForCSV() {
        return String.format("%s,%.2f,%.2f", toString(), getEstimatedQReported(), getEmpiricalQuality() - getEstimatedQReported());
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // increment methods
    //
    //---------------------------------------------------------------------------------------------------------------

    public final double getNumObservations() {
        return numObservations;
    }

    public final synchronized void setNumObservations(final double numObservations) {
        if ( numObservations < 0 ) throw new IllegalArgumentException("numObservations < 0");
        this.numObservations = numObservations;
        empiricalQuality = UNINITIALIZED;
    }

    public final double getNumMismatches() {
        return numMismatches;
    }

    @Requires({"numMismatches >= 0"})
    public final synchronized void setNumMismatches(final double numMismatches) {
        if ( numMismatches < 0 ) throw new IllegalArgumentException("numMismatches < 0");
        this.numMismatches = numMismatches;
        empiricalQuality = UNINITIALIZED;
    }

    @Requires({"by >= 0"})
    public final synchronized void incrementNumObservations(final double by) {
        numObservations += by;
        empiricalQuality = UNINITIALIZED;
    }

    @Requires({"by >= 0"})
    public final synchronized void incrementNumMismatches(final double by) {
        numMismatches += by;
        empiricalQuality = UNINITIALIZED;
    }

    @Requires({"incObservations >= 0", "incMismatches >= 0"})
    @Ensures({"numObservations == old(numObservations) + incObservations", "numMismatches == old(numMismatches) + incMismatches"})
    public final synchronized void increment(final double incObservations, final double incMismatches) {
        numObservations += incObservations;
        numMismatches += incMismatches;
        empiricalQuality = UNINITIALIZED;
    }

    @Ensures({"numObservations == old(numObservations) + 1", "numMismatches >= old(numMismatches)"})
    public final synchronized void increment(final boolean isError) {
        increment(1, isError ? 1 : 0.0);
    }

    // -------------------------------------------------------------------------------------
    //
    // Private implementation helper functions
    //
    // -------------------------------------------------------------------------------------

    /**
     * Calculate and cache the empirical quality score from mismatches and observations (expensive operation)
     */
    @Requires("empiricalQuality == UNINITIALIZED")
    @Ensures("empiricalQuality != UNINITIALIZED")
    private synchronized void calcEmpiricalQuality() {

        // TODO -- add code for Bayesian estimate of Qemp here

        final double empiricalQual = -10 * Math.log10(getEmpiricalErrorRate());
        empiricalQuality = Math.min(empiricalQual, (double) QualityUtils.MAX_RECALIBRATED_Q_SCORE);
    }

    /**
     * calculate the expected number of errors given the estimated Q reported and the number of observations
     * in this datum.
     *
     * @return a positive (potentially fractional) estimate of the number of errors
     */
    @Ensures("result >= 0.0")
    private double calcExpectedErrors() {
        return getNumObservations() * QualityUtils.qualToErrorProb(estimatedQReported);
    }

    static final boolean DEBUG = false;
    static final double RESOLUTION_BINS_PER_QUAL = 1.0;

    static public double bayesianEstimateOfEmpiricalQuality(final double nObservations, final double nErrors, final double QReported) {

        final int numBins = (QualityUtils.MAX_REASONABLE_Q_SCORE + 1) * (int)RESOLUTION_BINS_PER_QUAL;

        final double[] log10Posteriors = new double[numBins];

        for ( int bin = 0; bin < numBins; bin++ ) {

            final double QEmpOfBin = bin / RESOLUTION_BINS_PER_QUAL;

            log10Posteriors[bin] = log10QempPrior(QEmpOfBin, QReported) + log10Likelihood(QEmpOfBin, nObservations, nErrors);

            if ( DEBUG )
                System.out.println(String.format("bin = %d, Qreported = %f, nObservations = %f, nErrors = %f, posteriors = %f", bin, QReported, nObservations, nErrors, log10Posteriors[bin]));
        }

        if ( DEBUG )
            System.out.println(String.format("Qreported = %f, nObservations = %f, nErrors = %f", QReported, nObservations, nErrors));

        final double[] normalizedPosteriors = MathUtils.normalizeFromLog10(log10Posteriors);
        final int MLEbin = MathUtils.maxElementIndex(normalizedPosteriors);

        final double Qemp = MLEbin / RESOLUTION_BINS_PER_QUAL;
        return Qemp;
    }

    static final double[] log10QempPriorCache = new double[QualityUtils.MAX_GATK_USABLE_Q_SCORE + 1];
    static {
        // f(x) = a + b*exp(-((x - c)^2 / (2*d^2)))
        // Note that b is the height of the curve's peak, c is the position of the center of the peak, and d controls the width of the "bell".
        final double GF_a = 0.0;
        final double GF_b = 0.9;
        final double GF_c = 0.0;
        final double GF_d = 0.5;

        final GaussianFunction gaussian = new GaussianFunction(GF_a, GF_b, GF_c, GF_d);
        for ( int i = 0; i <= QualityUtils.MAX_GATK_USABLE_Q_SCORE; i++ )
            log10QempPriorCache[i] = Math.log10(gaussian.value((double) i));
    }

    static public double log10QempPrior(final double Qempirical, final double Qreported) {
        final int difference = Math.min(Math.abs((int) (Qempirical - Qreported)), QualityUtils.MAX_GATK_USABLE_Q_SCORE);
        if ( DEBUG )
            System.out.println(String.format("Qemp = %f, log10Priors = %f", Qempirical, log10QempPriorCache[difference]));
        return log10QempPriorCache[difference];
    }

    static public double log10Likelihood(final double Qempirical, final double nObservations, final double nErrors) {
        // this is just a straight binomial PDF
        double log10Prob = MathUtils.log10BinomialProbability((int)nObservations, (int)nErrors, QualityUtils.qualToErrorProbLog10((byte)(int)Qempirical));
        if ( log10Prob == Double.NEGATIVE_INFINITY )
            log10Prob = -Double.MAX_VALUE;

        if ( DEBUG )
            System.out.println(String.format("Qemp = %f, log10Likelihood = %f", Qempirical, log10Prob));
        return log10Prob;
    }
}