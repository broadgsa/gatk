/*
 * Copyright (c) 2011 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.SAMRecord;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 9/26/11
 */

public class ContextCovariate implements Covariate {

    final int CONTEXT_SIZE = 8;
    String allN = "";

    // Initialize any member variables using the command-line arguments passed to the walkers
    public void initialize( final RecalibrationArgumentCollection RAC ) {
        for( int iii = 0; iii < CONTEXT_SIZE; iii++ ) {
            allN += "N";
        }
    }

    public void getValues(SAMRecord read, Comparable[] comparable) {
        byte[] bases = read.getReadBases();
        for(int i = 0; i < read.getReadLength(); i++) {
            comparable[i] = ( i-CONTEXT_SIZE < 0 ? allN : new String(Arrays.copyOfRange(bases,i-CONTEXT_SIZE,i)) );
        }
    }

    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public final Comparable getValue( final String str ) {
        return str;
    }

}
