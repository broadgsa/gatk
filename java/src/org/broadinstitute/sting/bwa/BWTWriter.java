package org.broadinstitute.sting.bwa;

import org.broadinstitute.sting.utils.StingException;

import java.io.*;
import java.nio.ByteOrder;

/**
 * Writes an in-memory BWT to an outputstream.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWTWriter {
    /**
     * Input stream from which to read BWT data.
     */
    private final OutputStream outputStream;

    /**
     * Create a new BWT reader.
     * @param inputFile File in which the BWT is stored.
     */
    public BWTWriter( File inputFile ) {
        try {
            this.outputStream = new BufferedOutputStream(new FileOutputStream(inputFile));
        }
        catch( FileNotFoundException ex ) {
            throw new StingException("Unable to open output file", ex);
        }
    }

    /**
     * Write a BWT to the output stream.
     * @param bwt Transform to be written to the output stream.
     */
    public void write( BWT bwt ) {
        IntPackedOutputStream intPackedOutputStream = new IntPackedOutputStream(outputStream, ByteOrder.LITTLE_ENDIAN);
        BasePackedOutputStream basePackedOutputStream = new BasePackedOutputStream<Integer>(Integer.class, outputStream, ByteOrder.LITTLE_ENDIAN);

        try {
            intPackedOutputStream.write(bwt.inverseSA0);
            intPackedOutputStream.write(bwt.counts.toArray(true));

            for( SequenceBlock block: bwt.sequenceBlocks ) {
                intPackedOutputStream.write(block.occurrences.toArray(false));
                basePackedOutputStream.write(block.sequence);
            }

            // The last block is the last set of counts in the structure.
            intPackedOutputStream.write(bwt.counts.toArray(false));
        }
        catch( IOException ex ) {
            throw new StingException("Unable to read BWT from input stream.", ex);
        }
    }

    /**
     * Close the input stream.
     */
    public void close() {
        try {
            outputStream.close();
        }
        catch( IOException ex ) {
            throw new StingException("Unable to close input file", ex);
        }
    }
}
