package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

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

    public String platform;

    public CycleCovariate() { // empty constructor is required to instantiate covariate in CovariateCounterWalker and TableRecalibrationWalker
        platform = null;
    }

    public CycleCovariate(String _platform) {
        platform = _platform;
    }

    public Comparable getValue(SAMRecord read, int offset, char[] refBases) {
        //BUGBUG: assumes Solexa platform
        Integer cycle = offset;
        if( read.getReadNegativeStrandFlag() ) {
            cycle = read.getReadLength() - (offset + 1);
        }
        return cycle;
    }
    
    public Comparable getValue(String str) {
        return Integer.parseInt( str );
    }

    public int estimatedNumberOfBins() {
        return 100;
    }

    public String toString() {
        return "Cycle";
    }
}