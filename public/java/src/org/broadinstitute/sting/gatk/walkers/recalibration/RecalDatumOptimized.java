package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.List;

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

public class RecalDatumOptimized {

    protected long numObservations; // number of bases seen in total
    protected long numMismatches; // number of bases seen that didn't match the reference

    //---------------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    //---------------------------------------------------------------------------------------------------------------

    public RecalDatumOptimized() {
        numObservations = 0L;
        numMismatches = 0L;
    }

    public RecalDatumOptimized( final long _numObservations, final long _numMismatches) {
        numObservations = _numObservations;
        numMismatches = _numMismatches;
    }

    public RecalDatumOptimized( final RecalDatumOptimized copy ) {
        this.numObservations = copy.numObservations;
        this.numMismatches = copy.numMismatches;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // increment methods
    //
    //---------------------------------------------------------------------------------------------------------------

    public synchronized final void increment( final long incObservations, final long incMismatches ) {
        numObservations += incObservations;
        numMismatches += incMismatches;
    }

    public synchronized final void increment( final RecalDatumOptimized other ) {
        increment( other.numObservations, other.numMismatches );
    }

    public synchronized final void increment( final List<RecalDatumOptimized> data ) {
        for ( RecalDatumOptimized other : data ) {
            this.increment( other );
        }
    }

    public synchronized final void incrementBaseCounts( final byte curBase, final byte refBase ) {
        increment( 1, BaseUtils.simpleBaseToBaseIndex(curBase) == BaseUtils.simpleBaseToBaseIndex(refBase) ? 0 : 1 ); // increment takes num observations, then num mismatches
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // methods to derive empirical quality score
    //
    //---------------------------------------------------------------------------------------------------------------

    public final double empiricalQualDouble( final int smoothing, final double maxQual ) {
        final double doubleMismatches = (double) ( numMismatches + smoothing );
        final double doubleObservations = (double) ( numObservations + smoothing );
        double empiricalQual = -10 * Math.log10(doubleMismatches / doubleObservations);
        if (empiricalQual > maxQual) { empiricalQual = maxQual; }
        return empiricalQual;
    }
    public final double empiricalQualDouble() { return empiricalQualDouble( 0, QualityUtils.MAX_REASONABLE_Q_SCORE ); } // 'default' behavior is to use smoothing value of zero

    public final byte empiricalQualByte( final int smoothing ) {
        final double doubleMismatches = (double) ( numMismatches + smoothing );
        final double doubleObservations = (double) ( numObservations + smoothing );
        return QualityUtils.probToQual( 1.0 - doubleMismatches / doubleObservations ); // This is capped at Q40
    }
    public final byte empiricalQualByte() { return empiricalQualByte( 0 ); } // 'default' behavior is to use smoothing value of zero

    //---------------------------------------------------------------------------------------------------------------
    //
    // misc. methods
    //
    //---------------------------------------------------------------------------------------------------------------

    public final long getNumObservations() {
        return numObservations;
    }

    public final long getNumMismatches() {
        return numMismatches;
    }

    public final String outputToCSV( ) {
        return String.format( "%d,%d,%d", numObservations, numMismatches, (int)empiricalQualByte() );
    }
    public final String outputToCSV( final int smoothing ) {
        return String.format( "%d,%d,%d", numObservations, numMismatches, (int)empiricalQualByte(smoothing) );
    }
}
