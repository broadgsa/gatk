package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.QualityUtils;

/*
 * Copyright (c) 2010 The Broad Institute
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

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Jan 6, 2010
 *
 * An individual piece of recalibration data. Optimized for CountCovariates. Extras added to make TableRecalibration fast have been removed.
 * Each bin counts up the number of observations and the number of reference mismatches seen for that combination of covariates.
 */

public class Datum {

    long numObservations;                                                                                     // number of bases seen in total
    long numMismatches;                                                                                       // number of bases seen that didn't match the reference

    private static final int SMOOTHING_CONSTANT = 1;                                                                    // used when calculating empirical qualities to avoid division by zero


    //---------------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    //---------------------------------------------------------------------------------------------------------------

    public Datum() {
        numObservations = 0L;
        numMismatches = 0L;
    }

    public Datum(long numObservations, long numMismatches) {
        this.numObservations = numObservations;
        this.numMismatches = numMismatches;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // increment methods
    //
    //---------------------------------------------------------------------------------------------------------------

    synchronized void increment(final long incObservations, final long incMismatches) {
        numObservations += incObservations;
        numMismatches += incMismatches;
    }

    synchronized void increment(final boolean isError) {
        numObservations++;
        numMismatches += isError ? 1:0;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // methods to derive empirical quality score
    //
    //---------------------------------------------------------------------------------------------------------------

    double empiricalQualDouble() {
        final double doubleMismatches = (double) (numMismatches + SMOOTHING_CONSTANT);
        final double doubleObservations = (double) (numObservations + SMOOTHING_CONSTANT);
        final double empiricalQual = -10 * Math.log10(doubleMismatches / doubleObservations);
        return Math.min(empiricalQual, (double) QualityUtils.MAX_RECALIBRATED_Q_SCORE);
    }

    byte empiricalQualByte() {
        final double doubleMismatches = (double) (numMismatches);
        final double doubleObservations = (double) (numObservations);
        return QualityUtils.probToQual(1.0 - doubleMismatches / doubleObservations);                                    // This is capped at Q40
    }

    @Override
    public String toString() {
        return String.format("%d,%d,%d", numObservations, numMismatches, (int) empiricalQualByte());
    }

    @Override
    public boolean equals(Object o) {
        if (!(o instanceof Datum))
            return false;
        Datum other = (Datum) o;
        return numMismatches == other.numMismatches && numObservations == other.numObservations;
    }
}
