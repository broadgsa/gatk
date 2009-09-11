package org.broadinstitute.sting.bwa;

/**
 * Represents the Burrows-Wheeler Transform of a reference sequence.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWT {
    /**
     * Write an occurrence table after every SEQUENCE_BLOCK_SIZE bases.
     * For this implementation to behave correctly, SEQUENCE_BLOCK_SIZE % 8 == 0
     */
    public static final int SEQUENCE_BLOCK_SIZE = 128;    

    public final int inverseSA0;
    public final Counts counts;
    public final SequenceBlock[] sequenceBlocks;

    /**
     * Creates a new BWT with the given inverse SA, counts, and sequence (in ASCII).
     * @param inverseSA0 Inverse SA entry for the first element.  Will be missing from the BWT sequence.
     * @param counts Cumulative count of bases, in A,C,G,T order.
     * @param sequenceBlocks The full BWT sequence, sans the '$'.
     */
    public BWT( int inverseSA0, Counts counts, SequenceBlock[] sequenceBlocks ) {
        this.inverseSA0 = inverseSA0;
        this.counts = counts;
        this.sequenceBlocks = sequenceBlocks;
    }

    /**
     * Creates a new BWT with the given inverse SA, occurrences, and sequence (in ASCII).
     * @param inverseSA0 Inverse SA entry for the first element.  Will be missing from the BWT sequence.
     * @param counts Count of bases, in A,C,G,T order.
     * @param sequence The full BWT sequence, sans the '$'.
     */
    public BWT( int inverseSA0, Counts counts, byte[] sequence ) {
        this(inverseSA0,counts,generateSequenceBlocks(sequence));
    }

    /**
     * Extract the full sequence from the list of block.
     * @return The full BWT string as a byte array.
     */
    public byte[] getSequence() {
        byte[] sequence = new byte[counts.getTotal()];
        for( SequenceBlock block: sequenceBlocks )
            System.arraycopy(block.sequence,0,sequence,block.sequenceStart,block.sequenceLength);
        return sequence;
    }

    /**
     * The number of bases in the BWT as a whole.
     * @return Number of bases.
     */
    public int length() {
        return counts.getTotal();
    }

    /**
     * Create a set of sequence blocks from one long sequence.
     * @param sequence Sequence from which to derive blocks.
     * @return Array of sequence blocks containing data from the sequence.
     */
    private static SequenceBlock[] generateSequenceBlocks( byte[] sequence ) {
        Counts occurrences = new Counts();

        int numSequenceBlocks = PackUtils.numberOfPartitions(sequence.length,SEQUENCE_BLOCK_SIZE);
        SequenceBlock[] sequenceBlocks = new SequenceBlock[numSequenceBlocks];

        for( int block = 0; block < numSequenceBlocks; block++ ) {
            int blockStart = block*SEQUENCE_BLOCK_SIZE;
            int blockLength = Math.min(SEQUENCE_BLOCK_SIZE, sequence.length-blockStart);
            byte[] subsequence = new byte[blockLength];

            System.arraycopy(sequence,blockStart,subsequence,0,blockLength);

            sequenceBlocks[block] = new SequenceBlock(blockStart,blockLength,occurrences.clone(),subsequence);

            for( byte base: subsequence )
                occurrences.increment(Base.fromASCII(base));
        }

        return sequenceBlocks;
    }
}

/**
 * Models a block of bases within the BWT.
 */
class SequenceBlock {
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
    SequenceBlock( int sequenceStart, int sequenceLength, Counts occurrences, byte[] sequence ) {
        this.sequenceStart = sequenceStart;
        this.sequenceLength = sequenceLength;
        this.occurrences = occurrences;
        this.sequence = sequence;
    }
}
