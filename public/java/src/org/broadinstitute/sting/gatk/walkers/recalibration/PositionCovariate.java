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
 * Date: Nov 18, 2009
 *
 * The Position covariate. It is a lot like the Cycle covariate except it always returns the offset regardless of which platform the read came from.
 * This is the Solexa definition of machine cycle and the covariate that was always being used in the original version of the recalibrator.
 */

public class PositionCovariate implements ExperimentalCovariate {

    // Initialize any member variables using the command-line arguments passed to the walkers
    public void initialize( final RecalibrationArgumentCollection RAC ) {
    }

    // Used to pick out the covariate's value from attributes of the read
    public final Comparable getValue( final SAMRecord read, final int offset ) {
        int cycle = offset;
        if( read.getReadNegativeStrandFlag() ) {
            cycle = read.getReadLength() - (offset + 1);
        }
        return cycle;
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