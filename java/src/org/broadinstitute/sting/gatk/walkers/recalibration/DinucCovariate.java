package org.broadinstitute.sting.gatk.walkers.recalibration;

import java.util.HashMap;

import org.broadinstitute.sting.utils.BaseUtils;
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
 * Date: Nov 3, 2009
 *
 * The Dinucleotide covariate. This base and the one that came before it in the read, remembering to swap directions if negative strand read.
 * This covariate assumes that the bases have been swapped to their complement base counterpart if this is a negative strand read.
 * This assumption is made to speed up the code.
 */

public class DinucCovariate implements StandardCovariate {

    private static final byte NO_CALL = (byte)'N';
    private static final Dinuc NO_DINUC = new Dinuc(NO_CALL, NO_CALL);

    HashMap<Integer, Dinuc> dinucHashMap;

    // Initialize any member variables using the command-line arguments passed to the walkers
    public void initialize( final RecalibrationArgumentCollection RAC ) {
        final byte[] BASES = { (byte)'A', (byte)'C', (byte)'G', (byte)'T' };
        dinucHashMap = new HashMap<Integer, Dinuc>();
        for( byte byte1 : BASES ) {
            for( byte byte2: BASES ) {
                dinucHashMap.put( Dinuc.hashBytes(byte1, byte2), new Dinuc(byte1, byte2) ); // This might seem silly, but Strings are too slow
            }
        }
        // Add the "no dinuc" entry too
        dinucHashMap.put( Dinuc.hashBytes(NO_CALL, NO_CALL), NO_DINUC );
    }

    // Used to pick out the covariate's value from attributes of the read
    public final Comparable getValue( final SAMRecord read, final int offset ) {

        byte base;
        byte prevBase;
        final byte[] bases = read.getReadBases();
        // If this is a negative strand read then we need to reverse the direction for our previous base
        if( read.getReadNegativeStrandFlag() ) {
            // No dinuc at the beginning of the read
            if( offset == bases.length-1 ) {
                return NO_DINUC;
            }
            base = (byte)BaseUtils.simpleComplement( (char)(bases[offset]) );
            // Note: We are using the previous base in the read, not the previous base in the reference. This is done in part to be consistent with unmapped reads.
            prevBase = (byte)BaseUtils.simpleComplement( (char)(bases[offset + 1]) );
        } else {
            // No dinuc at the beginning of the read
            if( offset == 0 ) {
                return NO_DINUC;
            }
            base = bases[offset];
            // Note: We are using the previous base in the read, not the previous base in the reference. This is done in part to be consistent with unmapped reads.
            prevBase = bases[offset - 1];
        }

        // Make sure the previous base is good
        if( !BaseUtils.isRegularBase( prevBase ) ) {
            return NO_DINUC;
        }
        
        return dinucHashMap.get( Dinuc.hashBytes( prevBase, base ) );
    }

    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public final Comparable getValue( final String str ) {
        final Dinuc returnDinuc = dinucHashMap.get( Dinuc.hashBytes( (byte)str.charAt(0), (byte)str.charAt(1) ) );
        if( returnDinuc.compareTo(NO_DINUC) == 0 ) {
            return null;
        }
        return returnDinuc;

    }

    // Used to estimate the amount space required for the full data HashMap
    public final int estimatedNumberOfBins() {
        return 16;
    }
}