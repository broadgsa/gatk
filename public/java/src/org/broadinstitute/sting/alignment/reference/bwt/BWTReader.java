package org.broadinstitute.sting.alignment.reference.bwt;

import org.broadinstitute.sting.alignment.reference.packing.BasePackedInputStream;
import org.broadinstitute.sting.alignment.reference.packing.PackUtils;
import org.broadinstitute.sting.alignment.reference.packing.UnsignedIntPackedInputStream;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
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
            throw new ReviewedStingException("Unable to open input file", ex);
        }
    }

    /**
     * Read a BWT from the input stream.
     * @return The BWT stored in the input stream.
     */
    public BWT read() {
        UnsignedIntPackedInputStream uintPackedInputStream = new UnsignedIntPackedInputStream(inputStream, ByteOrder.LITTLE_ENDIAN);
        BasePackedInputStream basePackedInputStream = new BasePackedInputStream<Integer>(Integer.class, inputStream, ByteOrder.LITTLE_ENDIAN);

        long inverseSA0;
        long[] count;
        SequenceBlock[] sequenceBlocks;

        try {
            inverseSA0 = uintPackedInputStream.read();
            count = new long[PackUtils.ALPHABET_SIZE];
            uintPackedInputStream.read(count);

            long bwtSize = count[PackUtils.ALPHABET_SIZE-1];
            sequenceBlocks = new SequenceBlock[PackUtils.numberOfPartitions(bwtSize,BWT.SEQUENCE_BLOCK_SIZE)];
            
            for( int block = 0; block < sequenceBlocks.length; block++ ) {
                int sequenceStart = block* BWT.SEQUENCE_BLOCK_SIZE;
                int sequenceLength = (int)Math.min(BWT.SEQUENCE_BLOCK_SIZE,bwtSize-sequenceStart);

                long[] occurrences = new long[PackUtils.ALPHABET_SIZE];
                byte[] bwt = new byte[sequenceLength];

                uintPackedInputStream.read(occurrences);
                basePackedInputStream.read(bwt);

                sequenceBlocks[block] = new SequenceBlock(sequenceStart,sequenceLength,new Counts(occurrences,false),bwt);
            }
        }
        catch( IOException ex ) {
            throw new ReviewedStingException("Unable to read BWT from input stream.", ex);
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
            throw new ReviewedStingException("Unable to close input file", ex);
        }
    }
}
