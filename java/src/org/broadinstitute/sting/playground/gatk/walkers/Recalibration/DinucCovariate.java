package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.BaseUtils;

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

    public DinucCovariate() { // empty constructor is required to instantiate covariate in CovariateCounterWalker and TableRecalibrationWalker
    }
    
    public Comparable getValue(SAMRecord read, int offset, char[] refBases) {
        byte[] bases = read.getReadBases();
        char base = (char)bases[offset];
        char prevBase; // value set below
        // If it is a negative strand read then we need to get the complement base and reverse the direction for our previous base
        if( read.getReadNegativeStrandFlag() ) {
            base = BaseUtils.simpleComplement(base);
            if ( offset + 1 > bases.length - 1 ) { return Covariate.COVARIATE_ERROR; } // no prevBase because at the beginning of the read
            prevBase = BaseUtils.simpleComplement( (char)bases[offset + 1] );
        } else {
            if ( offset - 1 < 0 ) { return Covariate.COVARIATE_ERROR; } // no prevBase because at the beginning of the read
            prevBase = (char)bases[offset - 1];
        }
        // Check if bad base, probably an 'N'
        if( ( BaseUtils.simpleBaseToBaseIndex(prevBase) == -1 ) || ( BaseUtils.simpleBaseToBaseIndex(base) == -1 ) ) {
            return Covariate.COVARIATE_NULL; // CovariateCounterWalker will recognize that null means skip this particular location in the read
        } else {
            return String.format("%c%c", prevBase, base);
        }
    }
    
    public Comparable getValue(String str) {
    	return str;
    }

    public int estimatedNumberOfBins() {
        return 16;
    }

    public String toString() {
        return "Dinucleotide";
    }
}