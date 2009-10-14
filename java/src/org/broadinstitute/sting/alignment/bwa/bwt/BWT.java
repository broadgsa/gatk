package org.broadinstitute.sting.alignment.bwa.bwt;

import org.broadinstitute.sting.alignment.bwa.packing.PackUtils;
import org.broadinstitute.sting.utils.StingException;

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

    /**
     * The inverse SA, used as a placeholder for determining where the special EOL character sits.
     */
    protected final int inverseSA0;

    /**
     * Cumulative counts for the entire BWT.
     */
    protected final Counts counts;

    /**
     * The individual sequence blocks, modelling how they appear on disk.
     */
    protected final SequenceBlock[] sequenceBlocks;

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
     * Get the total counts of bases lexicographically smaller than the given base, for Ferragina and Manzini's search.
     * @param base The base.
     * @return Total counts for all bases lexicographically smaller than this base.
     */
    public int counts(byte base) {
        return counts.getCumulative(base);
    }

    /**
     * Get the total counts of bases lexicographically smaller than the given base, for Ferragina and Manzini's search.
     * @param base The base.
     * @param index The position to search within the BWT.
     * @return Total counts for all bases lexicographically smaller than this base.
     */
    public int occurrences(byte base,int index) {
        SequenceBlock block = getSequenceBlock(index);
        int position = getSequencePosition(index);
        int accumulator = block.occurrences.get(base);
        for(int i = 0; i <= position; i++) {
            if(base == block.sequence[i])
                accumulator++;
        }
        return accumulator;
    }

    /**
     * The number of bases in the BWT as a whole.
     * @return Number of bases.
     */
    public int length() {
        return counts.getTotal();
    }

    /**
     * Gets the base at a given position in the BWT.
     * @param index The index to use.
     * @return The base at that location.
     */
    protected byte getBase(int index) {
        if(index == inverseSA0)
            throw new StingException(String.format("Base at index %d does not have a text representation",index));

        SequenceBlock block = getSequenceBlock(index);
        int position = getSequencePosition(index);
        return block.sequence[position];
    }

    private SequenceBlock getSequenceBlock(int index) {
        // If the index is above the SA-1[0], remap it to the appropriate coordinate space.
        if(index > inverseSA0) index--;
        return sequenceBlocks[index/SEQUENCE_BLOCK_SIZE];
    }

    private int getSequencePosition(int index) {
        // If the index is above the SA-1[0], remap it to the appropriate coordinate space.
        if(index > inverseSA0) index--;
        return index%SEQUENCE_BLOCK_SIZE;
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
                occurrences.increment(base);
        }

        return sequenceBlocks;
    }
}
