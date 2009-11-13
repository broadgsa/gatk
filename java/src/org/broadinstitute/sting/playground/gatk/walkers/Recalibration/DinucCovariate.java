package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Pair;

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
 * The Dinucleotide covariate. This base and the one that came before it in the read, remembering to swap directions on cardinality if negative strand read.
 * This covariate will return null if there are bases that don't belong such as 'N' or 'X'.
 * This covariate will also return null if attempting to get previous base at the start of the read.
 */

public class DinucCovariate implements Covariate {

	private String returnString;
	
    public DinucCovariate() { // empty constructor is required to instantiate covariate in CovariateCounterWalker and TableRecalibrationWalker
    }
    
    public final Comparable getValue(final SAMRecord read, final int offset, final String readGroup, 
			 final byte[] quals, final char[] bases, final char refBase) {
    	
        char base = bases[offset];
        char prevBase = bases[offset - 1]; 
        // If it is a negative strand read then we need to reverse the direction for our previous base
        if( read.getReadNegativeStrandFlag() ) {
            prevBase = bases[offset + 1];
        } 
        final char[] charArray = {prevBase, base};
        returnString = String.valueOf(charArray);
        return returnString;
        //return String.format("%c%c", prevBase, base); // This return statement is too slow
    }
    
    public final Comparable getValue(final String str) {
    	return str;
    }

    public final int estimatedNumberOfBins() {
        return 16;
    }

    public String toString() {
        return "Dinucleotide";
    }
}