package org.broadinstitute.sting.alignment.reference.bwt;

import org.broadinstitute.sting.alignment.reference.packing.BasePackedOutputStream;
import org.broadinstitute.sting.alignment.reference.packing.UnsignedIntPackedOutputStream;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

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
     * Create a new BWT writer.
     * @param outputFile File in which the BWT is stored.
     */
    public BWTWriter( File outputFile ) {
        try {
            this.outputStream = new BufferedOutputStream(new FileOutputStream(outputFile));
        }
        catch( FileNotFoundException ex ) {
            throw new ReviewedStingException("Unable to open output file", ex);
        }
    }

    /**
     * Write a BWT to the output stream.
     * @param bwt Transform to be written to the output stream.
     */
    public void write( BWT bwt ) {
        UnsignedIntPackedOutputStream intPackedOutputStream = new UnsignedIntPackedOutputStream(outputStream, ByteOrder.LITTLE_ENDIAN);
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
            throw new ReviewedStingException("Unable to read BWT from input stream.", ex);
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
            throw new ReviewedStingException("Unable to close input file", ex);
        }
    }
}
