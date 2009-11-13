package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.QualityUtils;

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
 * Date: Nov 4, 2009
 *
 * The Minimum Neighborhood Quality Score covariate, originally described by Chris Hartl.
 * This covariate is the minimum base quality score in the read in a small window around the current base.
 */

public class MinimumNQSCovariate implements Covariate {

    private int windowReach; // how far in each direction from the current base to look

    public MinimumNQSCovariate() { // empty constructor is required to instantiate covariate in CovariateCounterWalker and TableRecalibrationWalker
        windowReach = 1; // window size = 3 was the best covariate according to Chris's analysis
    }

    public MinimumNQSCovariate(final int windowSize) {
        windowReach = windowSize / 2; // integer division
    }

    public final Comparable getValue(final SAMRecord read, final int offset, final String readGroup, 
			 final byte[] quals, final char[] bases, final char refBase) {
    	
    	// Loop over the list of base quality scores in the window and find the minimum
        int minQual = quals[offset];
        int minIndex = Math.max(offset - windowReach, 0);
        int maxIndex = Math.min(offset + windowReach, quals.length - 1);
        for ( int iii = minIndex; iii < maxIndex; iii++ ) {
            if( quals[iii] < minQual ) {
                minQual = quals[iii];
            }
        }
        return minQual;
    }
    
    public final Comparable getValue(final String str) {
        return (int)Integer.parseInt( str ); // cast to primitive int (as opposed to Integer Object) is required so that the return value from the two getValue methods hash to same thing
    }

    public final int estimatedNumberOfBins() {
        return 40;
    }

    public String toString() {
        return "Minimum Neighborhood Quality Score";
    }
}
