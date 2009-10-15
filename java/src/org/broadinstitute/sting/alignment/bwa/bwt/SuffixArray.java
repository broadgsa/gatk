package org.broadinstitute.sting.alignment.bwa.bwt;

import org.broadinstitute.sting.utils.StingException;

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
            throw new StingException("A BWT must be provided if the sequence interval is not 1");
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
}
