package org.broadinstitute.sting.oneoffprojects.walkers.IndelCountCovariates;

import net.sf.samtools.SAMRecord;
//g import org.broadinstitute.sting.gatk.walkers.recalibration.ExperimentalCovariate;
//g import org.broadinstitute.sting.gatk.walkers.recalibration.RecalibrationArgumentCollection;


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
 * Date: Dec 4, 2009
 *
 * The Homopolymer Run Covariate. This is the number of consecutive bases in the previous N that match the current base.
 * For example, if N = 10:
 * ATTGCCCCGTAAAAAAAAATA
 * 001001230001234567800
 */

public class HomopolymerCovariate implements ExperimentalCovariate {

    int numBack = 7;

    // Initialize any member variables using the command-line arguments passed to the walkers
    public void initialize( final RecalibrationArgumentCollection RAC ) {
        numBack = RAC.HOMOPOLYMER_NBACK;
    }

    // Used to pick out the covariate's value from attributes of the read
    public final Comparable getValue( final SAMRecord read, final int offset ) {

        // This block of code is for if you don't want to only count consecutive bases
        // ATTGCCCCGTAAAAAAAAATA
        // 001001231211234567819
        /*
        int numAgree = 0; // The number of bases that agree with you in the previous numBack bases of the read
        int startPos = 0;
        int stopPos = 0;
        byte[] bases = read.getReadBases();
        byte thisBase = bases[offset];
        if( !read.getReadNegativeStrandFlag() ) { // Forward direction
            startPos = Math.max(offset - numBack, 0);
            stopPos = Math.max(offset - 1, 0);
        } else { // Negative direction
            startPos = Math.min(offset + 2, bases.length);
            stopPos = Math.min(offset + numBack + 1, bases.length);
        }

        for( int iii = startPos; iii < stopPos; iii++ ) {
            if( bases[iii] == thisBase ) { numAgree++; }
        }
        */

        int numAgree = 0; // The number of consecutive bases that agree with you in the previous numBack bases of the read
        final byte[] bases = read.getReadBases();
        int iii = offset;
        if( !read.getReadNegativeStrandFlag() ) { // Forward direction
            while( iii <= bases.length-2 && bases[iii] == bases[iii+1] && numAgree < numBack ) {
                numAgree++;
                iii++;
            }
        } else { // Negative direction
            while( iii >= 1 && bases[iii] == bases[iii-1] && numAgree < numBack ) {
                numAgree++;
                iii--;
            }
        }

        return numAgree;
    }

    private void getContextHomopolymerLength(final byte[] refBytes, Comparable[] hrunArray) {
        // compute forward hrun length, example:
        // AGGTGACCCCCCTGAGAG
        // 001000012345000000
        int runCount = 0;
        hrunArray[0] = 0;
        int[] hforward = new int[hrunArray.length];
        int[] hreverse = new int[hrunArray.length];

        for (int i = 1; i < refBytes.length; i++) {
            if (refBytes[i] == refBytes[i-1])
                hforward[i] = hforward[i-1]+1;
            else
                hforward[i] = 0;
        }

        // do similar thing for reverse length, example:
        // AGGTGACCCCCCTGAGAG
        // 021000543210000000
        // and then accumulate with forward values.
        // Total:
        // AGGTGACCCCCCTGAGAG
        // 022000555555000000
        for (int i=refBytes.length-1; i > 0; i--) {
            if (refBytes[i-1] == refBytes[i])
                hreverse[i-1] += hreverse[i]+1;
        }

        for (int i = 1; i < refBytes.length; i++)
            hrunArray[i] = hforward[i]+hreverse[i];
    }

    public void getValues(SAMRecord read, Comparable[] comparable) {

    //    getContextHomopolymerLength(read.getReadBases(), comparable);
        for(int iii = 0; iii < read.getReadLength(); iii++) {
           comparable[iii] = getValue(read, iii); // BUGBUG: this can be optimized
        }
    }

    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public final Comparable getValue( final String str ) {
        return Integer.parseInt( str );
    }

}