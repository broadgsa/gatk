package org.broadinstitute.sting.alignment.bwa.bwt;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.alignment.bwa.packing.IntPackedInputStream;
import org.broadinstitute.sting.alignment.bwa.packing.BasePackedInputStream;
import org.broadinstitute.sting.alignment.bwa.packing.PackUtils;

import java.io.*;
import java.nio.ByteOrder;
/**
 * Reads a BWT from a given file.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWTReader {
    /**
     * Input stream from which to read BWT data.
     */
    private FileInputStream inputStream;

    /**
     * Create a new BWT reader.
     * @param inputFile File in which the BWT is stored.
     */
    public BWTReader( File inputFile ) {
        try {
            this.inputStream = new FileInputStream(inputFile);
        }
        catch( FileNotFoundException ex ) {
            throw new StingException("Unable to open input file", ex);
        }
    }

    /**
     * Read a BWT from the input stream.
     * @return The BWT stored in the input stream.
     */
    public BWT read() {
        IntPackedInputStream intPackedInputStream = new IntPackedInputStream(inputStream, ByteOrder.LITTLE_ENDIAN);
        BasePackedInputStream basePackedInputStream = new BasePackedInputStream<Integer>(Integer.class, inputStream, ByteOrder.LITTLE_ENDIAN);

        int inverseSA0;
        int[] count;
        SequenceBlock[] sequenceBlocks;

        try {
            inverseSA0 = intPackedInputStream.read();
            count = new int[PackUtils.ALPHABET_SIZE];
            intPackedInputStream.read(count);

            int bwtSize = count[PackUtils.ALPHABET_SIZE-1];
            sequenceBlocks = new SequenceBlock[PackUtils.numberOfPartitions(bwtSize,BWT.SEQUENCE_BLOCK_SIZE)];
            
            for( int block = 0; block < sequenceBlocks.length; block++ ) {
                int sequenceStart = block* BWT.SEQUENCE_BLOCK_SIZE;
                int sequenceLength = Math.min(BWT.SEQUENCE_BLOCK_SIZE,bwtSize-sequenceStart);

                int[] occurrences = new int[PackUtils.ALPHABET_SIZE];
                byte[] bwt = new byte[sequenceLength];

                intPackedInputStream.read(occurrences);
                basePackedInputStream.read(bwt);

                sequenceBlocks[block] = new SequenceBlock(sequenceStart,sequenceLength,new Counts(occurrences,false),bwt);
            }
        }
        catch( IOException ex ) {
            throw new StingException("Unable to read BWT from input stream.", ex);
        }

        return new BWT(inverseSA0, new Counts(count,true), sequenceBlocks);
    }

    /**
     * Close the input stream.
     */
    public void close() {
        try {
            inputStream.close();
        }
        catch( IOException ex ) {
            throw new StingException("Unable to close input file", ex);
        }
    }
}
