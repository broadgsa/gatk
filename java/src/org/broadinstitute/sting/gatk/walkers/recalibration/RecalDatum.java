package org.broadinstitute.sting.gatk.walkers.recalibration;

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

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 3, 2009
 *
 * An individual piece of recalibration data. Each bin counts up the number of observations and the number of reference mismatches seen for that combination of covariates. 
 */

public class RecalDatum extends RecalDatumOptimized {

    private double estimatedQReported; // estimated reported quality score based on combined data's individual q-reporteds and number of observations
    private double empiricalQuality; // the empirical quality for datums that have been collapsed together (by read group and reported quality, for example)

    //---------------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    //---------------------------------------------------------------------------------------------------------------

    public RecalDatum() {
        numObservations = 0L;
        numMismatches = 0L;
        estimatedQReported = 0.0;
        empiricalQuality = 0.0;
    }

    public RecalDatum( final long _numObservations, final long _numMismatches, final double _estimatedQReported, final double _empiricalQuality ) {
        numObservations = _numObservations;
        numMismatches = _numMismatches;
        estimatedQReported = _estimatedQReported;
        empiricalQuality = _empiricalQuality;
    }

    public RecalDatum( final RecalDatum copy ) {
        this.numObservations = copy.numObservations;
        this.numMismatches = copy.numMismatches;
        this.estimatedQReported = copy.estimatedQReported;
        this.empiricalQuality = copy.empiricalQuality;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // increment methods
    //
    //---------------------------------------------------------------------------------------------------------------

    public final void combine( final RecalDatum other ) {
        final double sumErrors = this.calcExpectedErrors() + other.calcExpectedErrors();
        this.increment( other.numObservations, other.numMismatches );
        this.estimatedQReported = -10 * Math.log10(sumErrors / (double)this.numObservations);
        //if( this.estimatedQReported > QualityUtils.MAX_REASONABLE_Q_SCORE ) { this.estimatedQReported = QualityUtils.MAX_REASONABLE_Q_SCORE; }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // methods to derive empirical quality score
    //
    //---------------------------------------------------------------------------------------------------------------

    public final void calcCombinedEmpiricalQuality( final int smoothing, final int maxQual ) {
        this.empiricalQuality = empiricalQualDouble(smoothing, maxQual); // cache the value so we don't call log over and over again
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // misc. methods
    //
    //---------------------------------------------------------------------------------------------------------------

    public final double getEstimatedQReported() {
        return estimatedQReported;
    }

    public final double getEmpiricalQuality() {
        return empiricalQuality;
    }

    private double calcExpectedErrors() {
        return (double)this.numObservations * qualToErrorProb( estimatedQReported );
    }

    private double qualToErrorProb( final double qual ) {
        return Math.pow(10.0, qual / -10.0);
    }
}
