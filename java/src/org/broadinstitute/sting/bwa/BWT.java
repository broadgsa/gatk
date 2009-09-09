package org.broadinstitute.sting.bwa;

/**
 * Represents the Burrows-Wheeler Transform of a reference sequence.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWT {
    public final int inverseSA0;
    public final int[] occurrences;
    public final byte[] sequence;

    /**
     * Creates a new BWT with the given inverse SA, occurrences, and sequence (in ASCII).
     * @param inverseSA0 Inverse SA entry for the first element.  Will be missing from the BWT sequence.
     * @param occurrences Cumulative number of occurrences of A,C,G,T, in order.
     * @param sequence The full BWT sequence, sans the '$'.
     */
    public BWT( int inverseSA0, int[] occurrences, byte[] sequence ) {
        this.inverseSA0 = inverseSA0;
        this.occurrences = occurrences;
        this.sequence = sequence;
    }
}
