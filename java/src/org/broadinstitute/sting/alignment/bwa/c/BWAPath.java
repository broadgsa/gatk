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

package org.broadinstitute.sting.alignment.bwa.c;

/**
 * Models a BWA path.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWAPath {
    /**
     * Number of mismatches encountered along this path.
     */
    public final int numMismatches;

    /**
     * Number of gap opens encountered along this path.
     */
    public final int numGapOpens;

    /**
     * Number of gap extensions along this path.
     */
    public final int numGapExtensions;

    /**
     * Whether this alignment was found on the positive or negative strand.
     */
    public final boolean negativeStrand;

    /**
     * Starting coordinate in the BWT.
     */
    public final long k;

    /**
     * Ending coordinate in the BWT.
     */
    public final long l;

    /**
     * The score of this path.  
     */
    public final int score;

    /**
     * The number of best alignments seen along this path.
     */
    public final int bestCount;

    /**
     * The number of second best alignments seen along this path.
     */
    public final int secondBestCount;

    /**
     * Create a new path with the given attributes.
     * @param numMismatches Number of mismatches along path.
     * @param numGapOpens Number of gap opens along path.
     * @param numGapExtensions Number of gap extensions along path.
     * @param k Index to first coordinate within BWT.
     * @param l Index to last coordinate within BWT.
     * @param score Score of this alignment.  Not the mapping quality.
     */
    public BWAPath(int numMismatches, int numGapOpens, int numGapExtensions, boolean negativeStrand, long k, long l, int score, int bestCount, int secondBestCount) {
        this.numMismatches = numMismatches;
        this.numGapOpens = numGapOpens;
        this.numGapExtensions = numGapExtensions;
        this.negativeStrand = negativeStrand;
        this.k = k;
        this.l = l;
        this.score = score;
        this.bestCount = bestCount;
        this.secondBestCount = secondBestCount;
    }

}
