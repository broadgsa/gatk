package org.broadinstitute.sting.gatk.walkers.bqsr;

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

import org.broadinstitute.sting.utils.MathUtils;

import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 3, 2009
 *
 * An individual piece of recalibration data. Each bin counts up the number of observations and the number of reference mismatches seen for that combination of covariates.
 */

public class RecalDatum extends Datum {

    private double estimatedQReported;                                                                                  // estimated reported quality score based on combined data's individual q-reporteds and number of observations
    private double empiricalQuality;                                                                                    // the empirical quality for datums that have been collapsed together (by read group and reported quality, for example)


    //---------------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    //---------------------------------------------------------------------------------------------------------------

    public RecalDatum() {
        numObservations = 0L;
        numMismatches = 0L;
        estimatedQReported = 0.0;
        empiricalQuality = -1.0;
    }

    public RecalDatum(final long _numObservations, final long _numMismatches, final double _estimatedQReported, final double _empiricalQuality) {
        numObservations = _numObservations;
        numMismatches = _numMismatches;
        estimatedQReported = _estimatedQReported;
        empiricalQuality = _empiricalQuality;
    }

    public RecalDatum(final RecalDatum copy) {
        this.numObservations = copy.numObservations;
        this.numMismatches = copy.numMismatches;
        this.estimatedQReported = copy.estimatedQReported;
        this.empiricalQuality = copy.empiricalQuality;
    }

    public void combine(final RecalDatum other) {
        final double sumErrors = this.calcExpectedErrors() + other.calcExpectedErrors();
        this.increment(other.numObservations, other.numMismatches);
        this.estimatedQReported = -10 * Math.log10(sumErrors / this.numObservations);
        this.empiricalQuality = -1.0;                                                                                   // reset the empirical quality calculation so we never have a wrongly calculated empirical quality stored
    }

    public final void calcCombinedEmpiricalQuality() {
        this.empiricalQuality = empiricalQualDouble();                                                                  // cache the value so we don't call log over and over again
    }
    
    public final void calcEstimatedReportedQuality() {
        this.estimatedQReported = -10 * Math.log10(calcExpectedErrors() / numObservations);
    }

    public final double getEstimatedQReported() {
        return estimatedQReported;
    }

    public final double getEmpiricalQuality() {
        if (empiricalQuality < 0)
            calcCombinedEmpiricalQuality();
        return empiricalQuality;
    }

    /**
     * Makes a hard copy of the recal datum element
     *
     * @return a new recal datum object with the same contents of this datum.
     */
    public RecalDatum copy() {
        return new RecalDatum(numObservations, numMismatches, estimatedQReported, empiricalQuality);
    }

    @Override
    public String toString() {
        return String.format("%d,%d,%d", numObservations, numMismatches, (byte) Math.floor(getEmpiricalQuality()));
    }

    public String stringForCSV() {
        return String.format("%s,%d,%.2f", toString(), (byte) Math.floor(getEstimatedQReported()), getEmpiricalQuality() - getEstimatedQReported());
    }

    private double calcExpectedErrors() {
        return (double) this.numObservations * qualToErrorProb(estimatedQReported);
    }

    private double qualToErrorProb(final double qual) {
        return Math.pow(10.0, qual / -10.0);
    }

    public static RecalDatum createRandomRecalDatum(int maxObservations, int maxErrors) {
        Random random = new Random();
        int nObservations = random.nextInt(maxObservations);
        int nErrors = random.nextInt(maxErrors);
        Datum datum = new Datum(nObservations, nErrors);
        double empiricalQuality = datum.empiricalQualDouble();
        double estimatedQReported = empiricalQuality + ((10 * random.nextDouble()) - 5);                                // empirical quality +/- 5.
        return new RecalDatum(nObservations, nErrors, estimatedQReported, empiricalQuality);
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
}
