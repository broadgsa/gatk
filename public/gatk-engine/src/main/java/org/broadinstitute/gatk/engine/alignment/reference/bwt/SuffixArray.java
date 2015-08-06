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

package org.broadinstitute.gatk.engine.alignment.reference.bwt;

import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.Comparator;
import java.util.TreeSet;

/**
 * An in-memory representation of a suffix array.
 *
 * @author mhanna
 * @version 0.1
 */
public class SuffixArray {
    public final long inverseSA0;
    public final Counts occurrences;

    /**
     * The elements of the sequence actually stored in memory.
     */
    protected final long[] sequence;

    /**
     * How often are individual elements in the sequence actually stored
     * in memory, as opposed to being calculated on the fly?
     */
    protected final int sequenceInterval;

    /**
     * The BWT used to calculate missing portions of the sequence.
     */
    protected final BWT bwt;

    public SuffixArray(long inverseSA0, Counts occurrences, long[] sequence) {
        this(inverseSA0,occurrences,sequence,1,null);
    }

    /**
     * Creates a new sequence array with the given inverse SA, occurrences, and values.
     * @param inverseSA0 Inverse SA entry for the first element.
     * @param occurrences Cumulative number of occurrences of A,C,G,T, in order.
     * @param sequence The full suffix array.
     * @param sequenceInterval How frequently is the sequence interval stored.
     * @param bwt bwt used to infer the remaining entries in the BWT.
     */
    public SuffixArray(long inverseSA0, Counts occurrences, long[] sequence, int sequenceInterval, BWT bwt) {
        this.inverseSA0 = inverseSA0;
        this.occurrences = occurrences;
        this.sequence = sequence;
        this.sequenceInterval = sequenceInterval;
        this.bwt = bwt;

        if(sequenceInterval != 1 && bwt == null)
            throw new ReviewedGATKException("A BWT must be provided if the sequence interval is not 1");
    }

    /**
     * Retrieves the length of the sequence array.
     * @return Length of the suffix array.
     */
    public long length() {
        if( bwt != null )
            return bwt.length()+1;
        else
            return sequence.length;
    }

    /**
     * Get the suffix array value at a given sequence.
     * @param index Index at which to retrieve the suffix array vaule.
     * @return The suffix array value at that entry.
     */
    public long get(long index) {
        int iterations = 0;
        while(index%sequenceInterval != 0) {
            // The inverseSA0 ('$') doesn't have a usable ASCII representation; it must be treated as a special case.
            if(index == inverseSA0)
                index = 0;
            else {
                byte base = bwt.getBase(index);
                index = bwt.counts(base) + bwt.occurrences(base,index);
            }
            iterations++;
        }
        return (sequence[(int)(index/sequenceInterval)]+iterations) % length();
    }

    /**
     * Create a suffix array from a given reference sequence.
     * @param sequence The reference sequence to use when building the suffix array.
     * @return a constructed suffix array.
     */
    public static SuffixArray createFromReferenceSequence(byte[] sequence) {
        // The builder for the suffix array.  Use an integer in this case because
        // Java arrays can only hold an integer.
        TreeSet<Integer> suffixArrayBuilder = new TreeSet<Integer>(new SuffixArrayComparator(sequence));

        Counts occurrences = new Counts();
        for( byte base: sequence )
            occurrences.increment(base);

        // Build out the suffix array using a custom comparator.
        for( int i = 0; i <= sequence.length; i++ )
            suffixArrayBuilder.add(i);

        // Copy the suffix array into an array.
        long[] suffixArray = new long[suffixArrayBuilder.size()];
        int i = 0;
        for( Integer element: suffixArrayBuilder )
            suffixArray[i++] = element;

        // Find the first element in the inverse suffix array.
        long inverseSA0 = -1;
        for(i = 0; i < suffixArray.length; i++) {
            if(suffixArray[i] == 0)
                inverseSA0 = i;
        }
        if(inverseSA0 < 0)
            throw new ReviewedGATKException("Unable to find first inverse SA entry in generated suffix array.");

        return new SuffixArray(inverseSA0,occurrences,suffixArray);
    }    

    /**
     * Compares two suffix arrays of the given sequence.  Will return whichever string appears
     * first in lexicographic order.
     */
    private static class SuffixArrayComparator implements Comparator<Integer> {
        /**
         * The data source for all suffix arrays.
         */
        private final String sequence;

        /**
         * Create a new comparator.
         * @param sequence Reference sequence to use as basis for comparison.
         */
        public SuffixArrayComparator( byte[] sequence ) {
            // Processing the suffix array tends to be easier as a string.
            this.sequence = StringUtil.bytesToString(sequence);
        }

        /**
         * Compare the two given suffix arrays.  Criteria for comparison is the lexicographic order of
         * the two substrings sequence[lhs:], sequence[rhs:].
         * @param lhs Left-hand side of comparison.
         * @param rhs Right-hand side of comparison.
         * @return How the suffix arrays represented by lhs, rhs compare.
         */
        public int compare( Integer lhs, Integer rhs ) {
            String lhsSuffixArray = sequence.substring(lhs);
            String rhsSuffixArray = sequence.substring(rhs);
            return lhsSuffixArray.compareTo(rhsSuffixArray);
        }
    }

}
