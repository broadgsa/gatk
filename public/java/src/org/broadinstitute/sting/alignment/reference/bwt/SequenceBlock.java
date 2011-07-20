package org.broadinstitute.sting.alignment.reference.bwt;

/**
 * Models a block of bases within the BWT.
 */
public class SequenceBlock {
    /**
     * Start position of this sequence within the BWT.
     */
    public final int sequenceStart;

    /**
     * Length of this sequence within the BWT.
     */
    public final int sequenceLength;


    /**
     * Occurrences of each letter up to this sequence block.
     */
    public final Counts occurrences;

    /**
     * Sequence for this segment.
     */
    public final byte[] sequence;

    /**
     * Create a new block within this BWT.
     * @param sequenceStart Starting position of this sequence within the BWT.
     * @param sequenceLength Length of this sequence.
     * @param occurrences How many of each base has been seen before this sequence began.
     * @param sequence The actual sequence from the BWT.
     */
    public SequenceBlock( int sequenceStart, int sequenceLength, Counts occurrences, byte[] sequence ) {
        this.sequenceStart = sequenceStart;
        this.sequenceLength = sequenceLength;
        this.occurrences = occurrences;
        this.sequence = sequence;
    }
}