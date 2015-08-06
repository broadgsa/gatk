/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.alignment.bwa.java;

import org.broadinstitute.gatk.engine.alignment.reference.bwt.BWT;

import java.util.ArrayList;
import java.util.List;

/**
 * At any point along the given read, what is a good lower bound for the
 * total number of differences?
 *
 * @author mhanna
 * @version 0.1
 */
public class LowerBound {
    /**
     * Lower bound of the suffix array.
     */
    public final long loIndex;

    /**
     * Upper bound of the suffix array.
     */
    public final long hiIndex;

    /**
     * Width of the bwt from loIndex -> hiIndex, inclusive.
     */
    public final long width;

    /**
     * The lower bound at the given point.
     */
    public final int value;

    /**
     * Create a new lower bound with the given value.
     * @param loIndex The lower bound of the BWT.
     * @param hiIndex The upper bound of the BWT.
     * @param value Value for the lower bound at this site.
     */
    private LowerBound(long loIndex, long hiIndex, int value) {
        this.loIndex = loIndex;
        this.hiIndex = hiIndex;
        this.width = hiIndex - loIndex + 1;
        this.value = value;
    }

    /**
     * Create a non-optimal bound according to the algorithm specified in Figure 3 of the BWA paper.
     * @param bases Bases of the read to use when creating a new BWT.
     * @param bwt BWT to check against.
     * @return A list of lower bounds at every point in the reference.
     *
     */
    public static List<LowerBound> create(Byte[] bases, BWT bwt) {
        List<LowerBound> bounds = new ArrayList<LowerBound>();

        long loIndex = 0, hiIndex = bwt.length();
        int mismatches = 0;
        for( int i = bases.length-1; i >= 0; i-- ) {
            Byte base = bases[i];

            // Ignore non-ACGT bases.
            if( base != null ) {
                loIndex = bwt.counts(base) + bwt.occurrences(base,loIndex-1) + 1;
                hiIndex = bwt.counts(base) + bwt.occurrences(base,hiIndex);            
            }

            if( base == null || loIndex > hiIndex ) {
                loIndex = 0;
                hiIndex = bwt.length();
                mismatches++;
            }
            bounds.add(0,new LowerBound(loIndex,hiIndex,mismatches));
        }

        return bounds;
    }

    /**
     * Create a string representation of this bound.
     * @return String version of this bound.
     */
    public String toString() {
        return String.format("LowerBound: w = %d, value = %d",width,value);
    }
}
