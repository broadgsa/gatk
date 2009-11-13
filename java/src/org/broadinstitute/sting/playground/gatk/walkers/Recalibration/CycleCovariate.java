package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.utils.StingException;

import net.sf.samtools.SAMRecord;

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
 *  Yet to be implemented:
 *  - For 454 the cycle is the number of discontinuous nucleotides seen during the length of the read
 *     For example, for the read: AAACCCCGAAATTTTTTT
 *             the cycle would be 111222234445555555
 *  - For SOLiD the cycle is a more complicated mixture of ligation cycle and primer round 
 */

public class CycleCovariate implements Covariate {

	private String platform;

    public CycleCovariate() { // empty constructor is required to instantiate covariate in CovariateCounterWalker and TableRecalibrationWalker
        platform = "SLX";
    }

    public CycleCovariate(final String _platform) {
    	platform = _platform;
    }

    public final Comparable getValue(final SAMRecord read, final int offset, final String readGroup, 
    								 final byte[] quals, final char[] bases, final char refBase) {
        if( platform.equalsIgnoreCase( "SLX" ) ) {
	        int cycle = offset;
	        if( read.getReadNegativeStrandFlag() ) {
	            cycle = bases.length - (offset + 1);
	        }
	        return cycle;
        } else if( platform.equalsIgnoreCase( "454" ) ) {
        	return 0; //BUGBUG: not yet implemented
        } else if( platform.equalsIgnoreCase( "SOLID" ) ) {
        	return 0; //BUGBUG: not yet implemented
        } else {
        	throw new StingException( "Requested platform (" + platform + ") not supported in CycleCovariate." );
        }
                                                        
    }
    
    public final Comparable getValue(final String str) {
        return (int)Integer.parseInt( str );
    }

    public final int estimatedNumberOfBins() {
        return 100;
    }

    public String toString() {
        return "Cycle";
    }
}