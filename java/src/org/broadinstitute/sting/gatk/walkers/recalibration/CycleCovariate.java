package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.utils.Utils;

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
 * Date: Oct 30, 2009
 *
 * The Cycle covariate.
 *  For Solexa the cycle is simply the position in the read (counting backwards if it is a negative strand read)
 *  For 454 the cycle is the number of discontinuous nucleotides seen during the length of the read
 *     For example, for the read: AAACCCCGAAATTTTTTT
 *             the cycle would be 111222234445555555
 *  For SOLiD the cycle is a more complicated mixture of ligation cycle and primer round
 */

public class CycleCovariate implements Covariate {

    private static boolean warnedUserNoPlatform = false;

    public CycleCovariate() { // empty constructor is required to instantiate covariate in CovariateCounterWalker and TableRecalibrationWalker
    }

    public final Comparable getValue( final ReadHashDatum readDatum, final int offset ) {
        if( readDatum.platform.equalsIgnoreCase( "ILLUMINA" ) ) {
            int cycle = offset;
	        if( readDatum.isNegStrand ) {
	            cycle = readDatum.bases.length - (offset + 1);
	        }
	        return cycle;
        } else if( readDatum.platform.contains( "454" ) ) { // some bams have "LS454" and others have just "454"
            int cycle = 0;
            byte prevBase = readDatum.bases[0];
            for( int iii = 1; iii <= offset; iii++ ) {
                if(readDatum.bases[iii] != prevBase) { // this base doesn't match the previous one so it is a new cycle
                    cycle++;
                    prevBase = readDatum.bases[iii];
                }
            }
            return cycle;
        } else if( readDatum.platform.equalsIgnoreCase( "SOLID" ) ) {
            // the ligation cycle according to http://www3.appliedbiosystems.com/cms/groups/mcb_marketing/documents/generaldocuments/cms_057511.pdf
        	return offset / 5; // integer division
        } else { // platform is unrecognized so revert to Illumina definition of cycle but warn the user
        	if( !warnedUserNoPlatform ) {
                Utils.warnUser( "Platform (" + readDatum.platform + ") unrecognized. Reverting to Illumina definition of machine cycle." );
                warnedUserNoPlatform = true;
            }
            return PositionCovariate.revertToPositionAsCycle( readDatum, offset );
        }
                                                        
    }
    
    public final Comparable getValue( final String str ) {
        return (int)Integer.parseInt( str ); // cast to primitive int (as opposed to Integer Object) is required so that the return value from the two getValue methods hash to same thing
    }

    public final int estimatedNumberOfBins() {
        return 100;
    }

    public String toString() {
        return "Machine Cycle";
    }
}