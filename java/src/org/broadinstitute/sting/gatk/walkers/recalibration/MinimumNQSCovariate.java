package org.broadinstitute.sting.gatk.walkers.recalibration;

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
 * Date: Nov 4, 2009
 *
 * The Minimum Neighborhood Quality Score covariate, originally described by Chris Hartl.
 * This covariate is the minimum base quality score in the read in a small window around the current base.
 */

public class MinimumNQSCovariate implements ExperimentalCovariate {

    private int windowReach; // How far in each direction from the current base to look

    // Initialize any member variables using the command-line arguments passed to the walkers
    public void initialize( final RecalibrationArgumentCollection RAC ) {
        windowReach = RAC.WINDOW_SIZE / 2; // integer division
    }

    // Used to pick out the covariate's value from attributes of the read
    public final Comparable getValue( final SAMRecord read, final int offset ) {

        // Loop over the list of base quality scores in the window and find the minimum
        final byte[] quals = read.getBaseQualities();
        int minQual = quals[offset];
        final int minIndex = Math.max(offset - windowReach, 0);
        final int maxIndex = Math.min(offset + windowReach, quals.length - 1);
        for ( int iii = minIndex; iii < maxIndex; iii++ ) {
            if( quals[iii] < minQual ) {
                minQual = quals[iii];
            }
        }
        return minQual;
    }

    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public final Comparable getValue( final String str ) {
        return Integer.parseInt( str );
    }

    public void getValues(SAMRecord read, Comparable[] comparable) {
        for(int iii = 0; iii < read.getReadLength(); iii++) {
            comparable[iii] = getValue(read, iii); // BUGBUG: this can be optimized
        }
    }
}
