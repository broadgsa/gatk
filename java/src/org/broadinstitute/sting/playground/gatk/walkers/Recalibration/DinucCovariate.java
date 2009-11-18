package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

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
 * The Dinucleotide covariate. This base and the one that came before it in the read, remembering to swap directions if negative strand read.
 * This covariate assumes that the bases have been swapped to their complement base counterpart if this is a negative strand read.
 * This assumption is made to speed up the code.
 */

public class DinucCovariate implements Covariate {

    HashMap<Integer, Dinuc> dinucHashMap;

    public DinucCovariate() { // empty constructor is required to instantiate covariate in CovariateCounterWalker and TableRecalibrationWalker

        final byte[] BASES = { (byte)'A', (byte)'C', (byte)'G', (byte)'T' };
        dinucHashMap = new HashMap<Integer, Dinuc>();
        for(byte byte1 : BASES) {
            for(byte byte2: BASES) {
                dinucHashMap.put( Dinuc.hashBytes(byte1, byte2), new Dinuc(byte1, byte2) ); // This might seem silly, but Strings are too slow
            }
        }
    }
    
    public final Comparable getValue(final SAMRecord read, final int offset, final String readGroup, 
			 final byte[] quals, final byte[] bases) {
    	
        byte base = bases[offset];
        byte prevBase = bases[offset - 1];
        // If this is a negative strand read then we need to reverse the direction for our previous base
        if( read.getReadNegativeStrandFlag() ) {
            prevBase = bases[offset + 1];
        }
        //char[] charArray = {(char)prevBase, (char)base};
        //return new String( charArray ); // This is an expensive call
        return dinucHashMap.get( Dinuc.hashBytes( prevBase, base ) );
        //return String.format("%c%c", prevBase, base); // This return statement is too slow
    }
    
    public final Comparable getValue(final String str) {
        //return str;
        return dinucHashMap.get( Dinuc.hashBytes( (byte)str.charAt(0), (byte)str.charAt(1) ) );
    }

    public final int estimatedNumberOfBins() {
        return 16;
    }

    public String toString() {
        return "Dinucleotide";
    }
}