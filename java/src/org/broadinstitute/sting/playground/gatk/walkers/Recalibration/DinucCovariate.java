package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.BaseUtils;

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
 */

public class DinucCovariate implements Covariate {

    public static ArrayList<String> BASES;

    public DinucCovariate() { // empty constructor is required to instantiate covariate in CovariateCounterWalker and TableRecalibrationWalker
        BASES = new ArrayList<String>();
        BASES.add("A");
        BASES.add("G");
        BASES.add("C");
        BASES.add("T");
        BASES.add("a");
        BASES.add("g");
        BASES.add("c");
        BASES.add("t");
    }
    
    public Comparable getValue(SAMRecord read, int offset, char[] refBases) {
        byte[] bases = read.getReadBases();
        char base = (char)bases[offset];
        char prevBase = (char)bases[offset - 1];
        if( read.getReadNegativeStrandFlag() ) {
            base = BaseUtils.simpleComplement(base);
            prevBase = BaseUtils.simpleComplement( (char)bases[offset + 1] );
        }
        // Check if bad base, probably an 'N'
        if( !BASES.contains( String.format( "%c", prevBase ) ) || !BASES.contains( String.format( "%c", base) ) ) {
            return null; // CovariateCounterWalker will recognize that null means skip this particular location in the read
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