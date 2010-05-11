package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.SAMRecord;

/*
 * Copyright (c) 2010 The Broad Institute
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
 * Date: Jan 29, 2010
 *
 * The number of previous N bases (along the direction of the read) that contain G's and C's. The goal is to correct for dye slippage.
 * Only valid for Illumina reads. Otherwise return -1.
 */

public class GCContentCovariate implements ExperimentalCovariate {

    int numBack = 7;

    // Initialize any member variables using the command-line arguments passed to the walkers
    public void initialize( final RecalibrationArgumentCollection RAC ) {
        numBack = RAC.HOMOPOLYMER_NBACK;
    }

    // Used to pick out the covariate's value from attributes of the read
    public final Comparable getValue( final SAMRecord read, final int offset ) {

        // ATTGCCCCGTAAAAAAAGAGAA
        // 0000123456654321001122

        if( read.getReadGroup().getPlatform().equalsIgnoreCase( "ILLUMINA" ) || read.getReadGroup().getPlatform().equalsIgnoreCase( "SLX" ) ) {
            int numGC = 0;
            int startPos = 0;
            int stopPos = 0;
            final byte[] bases = read.getReadBases();
            if( !read.getReadNegativeStrandFlag() ) { // Forward direction
                startPos = Math.max(offset - numBack, 0);
                stopPos = Math.max(offset - 1, 0);
            } else { // Negative direction
                startPos = Math.min(offset + 2, bases.length);
                stopPos = Math.min(offset + numBack + 1, bases.length);
            }

            for( int iii = startPos; iii < stopPos; iii++ ) {
                if( bases[iii] == (byte)'G' || bases[iii] == (byte)'C' ) {
                    numGC++;
                }
            }

            return numGC;
        } else { // This effect is specific to the Illumina platform
            return -1;
        }
    }
    
    public void getValues(SAMRecord read, Comparable[] comparable) {
        for(int iii = 0; iii < read.getReadLength(); iii++) {
            comparable[iii] = getValue(read, iii); // BUGBUG: this can be optimized
        }
    }

    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public final Comparable getValue( final String str ) {
        return Integer.parseInt( str );
    }



}