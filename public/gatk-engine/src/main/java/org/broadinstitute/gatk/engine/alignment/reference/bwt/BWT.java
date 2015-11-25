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

import org.broadinstitute.gatk.engine.alignment.reference.packing.PackUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

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
    protected final long inverseSA0;

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
    public BWT( long inverseSA0, Counts counts, SequenceBlock[] sequenceBlocks ) {
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
    public BWT( long inverseSA0, Counts counts, byte[] sequence ) {
        this(inverseSA0,counts,generateSequenceBlocks(sequence));
    }

    /**
     * Extract the full sequence from the list of block.
     * @return The full BWT string as a byte array.
     */
    public byte[] getSequence() {
        byte[] sequence = new byte[(int)counts.getTotal()];
        for( SequenceBlock block: sequenceBlocks )
            System.arraycopy(block.sequence,0,sequence,block.sequenceStart,block.sequenceLength);
        return sequence;
    }

    /**
     * Get the total counts of bases lexicographically smaller than the given base, for Ferragina and Manzini's search.
     * @param base The base.
     * @return Total counts for all bases lexicographically smaller than this base.
     */
    public long counts(byte base) {
        return counts.getCumulative(base);
    }

    /**
     * Get the total counts of bases lexicographically smaller than the given base, for Ferragina and Manzini's search.
     * @param base The base.
     * @param index The position to search within the BWT.
     * @return Total counts for all bases lexicographically smaller than this base.
     */
    public long occurrences(byte base,long index) {
        SequenceBlock block = getSequenceBlock(index);
        int position = getSequencePosition(index);
        long accumulator = block.occurrences.get(base);
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
    public long length() {
        return counts.getTotal();
    }

    /**
     * Create a new BWT from the given reference sequence.
     * @param referenceSequence Sequence from which to derive the BWT.
     * @return reference sequence-derived BWT.
     */
    public static BWT createFromReferenceSequence(byte[] referenceSequence) {
        SuffixArray suffixArray = SuffixArray.createFromReferenceSequence(referenceSequence);

        byte[] bwt = new byte[(int)suffixArray.length()-1];
        int bwtIndex = 0;
        for(long suffixArrayIndex = 0; suffixArrayIndex < suffixArray.length(); suffixArrayIndex++) {
            if(suffixArray.get(suffixArrayIndex) == 0)
                continue;
            bwt[bwtIndex++] = referenceSequence[(int)suffixArray.get(suffixArrayIndex)-1];
        }

        return new BWT(suffixArray.inverseSA0,suffixArray.occurrences,bwt);
    }

    /**
     * Gets the base at a given position in the BWT.
     * @param index The index to use.
     * @return The base at that location.
     */
    protected byte getBase(long index) {
        if(index == inverseSA0)
            throw new ReviewedGATKException(String.format("Base at index %d does not have a text representation",index));

        SequenceBlock block = getSequenceBlock(index);
        int position = getSequencePosition(index);
        return block.sequence[position];
    }

    private SequenceBlock getSequenceBlock(long index) {
        // If the index is above the SA-1[0], remap it to the appropriate coordinate space.
        if(index > inverseSA0) index--;
        return sequenceBlocks[(int)(index/SEQUENCE_BLOCK_SIZE)];
    }

    private int getSequencePosition(long index) {
        // If the index is above the SA-1[0], remap it to the appropriate coordinate space.
        if(index > inverseSA0) index--;
        return (int)(index%SEQUENCE_BLOCK_SIZE);
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
