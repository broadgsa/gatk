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
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.Random;

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
    private long numObservations;

    /**
     * number of bases seen that didn't match the reference
     */
    private long numMismatches;

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
    public RecalDatum(final long _numObservations, final long _numMismatches, final byte reportedQuality) {
        if ( numObservations < 0 ) throw new IllegalArgumentException("numObservations < 0");
        if ( numMismatches < 0 ) throw new IllegalArgumentException("numMismatches < 0");
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

    @Override
    public String toString() {
        return String.format("%d,%d,%d", getNumObservations(), getNumMismatches(), (byte) Math.floor(getEmpiricalQuality()));
    }

    public String stringForCSV() {
        return String.format("%s,%d,%.2f", toString(), (byte) Math.floor(getEstimatedQReported()), getEmpiricalQuality() - getEstimatedQReported());
    }

    public static RecalDatum createRandomRecalDatum(int maxObservations, int maxErrors) {
        final Random random = new Random();
        final int nObservations = random.nextInt(maxObservations);
        final int nErrors = random.nextInt(maxErrors);
        final int qual = random.nextInt(QualityUtils.MAX_QUAL_SCORE);
        return new RecalDatum(nObservations, nErrors, (byte)qual);
    }

    /**
     * We don't compare the estimated quality reported because it may be different when read from
     * report tables.
     *
     * @param o the other recal datum
     * @return true if the two recal datums have the same number of observations, errors and empirical quality.
     */
    @Override
    public boolean equals(Object o) {
        if (!(o instanceof RecalDatum))
            return false;
        RecalDatum other = (RecalDatum) o;
        return super.equals(o) &&
               MathUtils.compareDoubles(this.empiricalQuality, other.empiricalQuality, 0.001) == 0;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // increment methods
    //
    //---------------------------------------------------------------------------------------------------------------

    public long getNumObservations() {
        return numObservations;
    }

    public synchronized void setNumObservations(final long numObservations) {
        if ( numObservations < 0 ) throw new IllegalArgumentException("numObservations < 0");
        this.numObservations = numObservations;
        empiricalQuality = UNINITIALIZED;
    }

    public long getNumMismatches() {
        return numMismatches;
    }

    @Requires({"numMismatches >= 0"})
    public synchronized void setNumMismatches(final long numMismatches) {
        if ( numMismatches < 0 ) throw new IllegalArgumentException("numMismatches < 0");
        this.numMismatches = numMismatches;
        empiricalQuality = UNINITIALIZED;
    }

    @Requires({"by >= 0"})
    public synchronized void incrementNumObservations(final long by) {
        numObservations += by;
        empiricalQuality = UNINITIALIZED;
    }

    @Requires({"by >= 0"})
    public synchronized void incrementNumMismatches(final long by) {
        numMismatches += by;
        empiricalQuality = UNINITIALIZED;
    }

    @Requires({"incObservations >= 0", "incMismatches >= 0"})
    @Ensures({"numObservations == old(numObservations) + incObservations", "numMismatches == old(numMismatches) + incMismatches"})
    public synchronized void increment(final long incObservations, final long incMismatches) {
        incrementNumObservations(incObservations);
        incrementNumMismatches(incMismatches);
    }

    @Ensures({"numObservations == old(numObservations) + 1", "numMismatches >= old(numMismatches)"})
    public synchronized void increment(final boolean isError) {
        incrementNumObservations(1);
        if ( isError )
            incrementNumMismatches(1);
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
    private synchronized final void calcEmpiricalQuality() {
        // cache the value so we don't call log over and over again
        final double doubleMismatches = (double) (numMismatches + SMOOTHING_CONSTANT);
        // smoothing is one error and one non-error observation, for example
        final double doubleObservations = (double) (numObservations + SMOOTHING_CONSTANT + SMOOTHING_CONSTANT);
        final double empiricalQual = -10 * Math.log10(doubleMismatches / doubleObservations);
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
        return (double) getNumObservations() * QualityUtils.qualToErrorProb(estimatedQReported);
    }
}
