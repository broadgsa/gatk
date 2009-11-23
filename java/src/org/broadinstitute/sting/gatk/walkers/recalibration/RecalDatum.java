package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.*;

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

public class RecalDatum {

    private long numObservations; // number of bases seen in total
    private long numMismatches; // number of bases seen that didn't match the reference


    //---------------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    //---------------------------------------------------------------------------------------------------------------
    public RecalDatum() {
        numObservations = 0L;
        numMismatches = 0L;
    }

    public RecalDatum( final long _numObservations, final long _numMismatches ) {
        numObservations = _numObservations;
        numMismatches = _numMismatches;
    }

    public RecalDatum( final RecalDatum copy ) {
        this.numObservations = copy.numObservations;
        this.numMismatches = copy.numMismatches;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // increment methods
    //
    //---------------------------------------------------------------------------------------------------------------


    public final void increment( final long incObservations, final long incMismatches ) {
        numObservations += incObservations;
        numMismatches += incMismatches;
    }

    public final void increment( final RecalDatum other ) {
        increment( other.numObservations, other.numMismatches );
    }

    public final void increment( final List<RecalDatum> data ) {
        for ( RecalDatum other : data ) {
            this.increment( other );
        }
    }

    public final void increment( final char curBase, final char refBase ) {
        increment( 1, BaseUtils.simpleBaseToBaseIndex(curBase) == BaseUtils.simpleBaseToBaseIndex(refBase) ? 0 : 1 ); // increment takes num observations, then num mismatches
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // methods to derive empirical quality score
    //
    //---------------------------------------------------------------------------------------------------------------

    public final double empiricalQualDouble( final int smoothing ) {
        double doubleMismatches = (double) ( numMismatches + smoothing );
        double doubleObservations = (double) ( numObservations + smoothing );
        double empiricalQual = -10 * Math.log10(doubleMismatches / doubleObservations);
        if (empiricalQual > QualityUtils.MAX_REASONABLE_Q_SCORE) { empiricalQual = QualityUtils.MAX_REASONABLE_Q_SCORE; }
        return empiricalQual;
    }
    public final double empiricalQualDouble() { return empiricalQualDouble( 0 ); } // 'default' behavior is to use smoothing value of zero


    public final byte empiricalQualByte( final int smoothing ) {
        double doubleMismatches = (double) ( numMismatches + smoothing );
        double doubleObservations = (double) ( numObservations + smoothing );
        return QualityUtils.probToQual( 1.0 - doubleMismatches / doubleObservations );
    }
    public final byte empiricalQualByte() { return empiricalQualByte( 0 ); } // 'default' behavior is to use smoothing value of zero

    //---------------------------------------------------------------------------------------------------------------
    //
    // misc. methods
    //
    //---------------------------------------------------------------------------------------------------------------

    public final String outputToCSV( ) {
        return String.format( "%d,%d,%d", numObservations, numMismatches, (int)empiricalQualByte() );
    }
    public final String outputToCSV( final int smoothing ) {
        return String.format( "%d,%d,%d", numObservations, numMismatches, (int)empiricalQualByte( smoothing ) );
    }

    public final Long getNumObservations() {
        return numObservations;
    }
    
    public String toString() {
    	return String.format( "RecalDatum: %d,%d,%d", numObservations, numMismatches, (int)empiricalQualByte() );
    }
}
