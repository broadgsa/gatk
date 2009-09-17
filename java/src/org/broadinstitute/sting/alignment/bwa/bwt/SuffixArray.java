package org.broadinstitute.sting.alignment.bwa.bwt;

/**
 * An in-memory representation of a suffix array.
 *
 * @author mhanna
 * @version 0.1
 */
public class SuffixArray {
    public final int inverseSA0;
    public final Counts occurrences;
    public final int[] sequence;

    /**
     * Creates a new sequence array with the given inverse SA, occurrences, and values.
     * @param inverseSA0 Inverse SA entry for the first element.
     * @param occurrences Cumulative number of occurrences of A,C,G,T, in order.
     * @param sequence The full suffix array.
     */
    public SuffixArray(int inverseSA0, Counts occurrences, int[] sequence) {
        this.inverseSA0 = inverseSA0;
        this.occurrences = occurrences;
        this.sequence = sequence;
    }

    /**
     * Retrieves the length of the sequence array.
     * @return Length of the suffix array.
     */
    public int length() {
        return sequence.length;
    }

}
