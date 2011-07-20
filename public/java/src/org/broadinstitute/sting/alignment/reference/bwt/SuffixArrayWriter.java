package org.broadinstitute.sting.alignment.reference.bwt;

import org.broadinstitute.sting.alignment.reference.packing.UnsignedIntPackedOutputStream;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.*;
import java.nio.ByteOrder;

/**
 * Javadoc goes here.
 *
 * @author mhanna
 * @version 0.1
 */
public class SuffixArrayWriter {
    /**
     * Input stream from which to read suffix array data.
     */
    private OutputStream outputStream;

    /**
     * Create a new suffix array reader.
     * @param outputFile File in which the suffix array is stored.
     */
    public SuffixArrayWriter( File outputFile ) {
        try {
            this.outputStream = new BufferedOutputStream(new FileOutputStream(outputFile));
        }
        catch( FileNotFoundException ex ) {
            throw new ReviewedStingException("Unable to open input file", ex);
        }
    }

    /**
     * Write a suffix array to the output stream.
     * @param suffixArray suffix array to write.
     */
    public void write(SuffixArray suffixArray) {
        UnsignedIntPackedOutputStream uintPackedOutputStream = new UnsignedIntPackedOutputStream(outputStream, ByteOrder.LITTLE_ENDIAN);

        try {
            uintPackedOutputStream.write(suffixArray.inverseSA0);
            uintPackedOutputStream.write(suffixArray.occurrences.toArray(true));
            // How frequently the suffix array entry is placed.
            uintPackedOutputStream.write(1);
            // Length of the suffix array.
            uintPackedOutputStream.write(suffixArray.length()-1);
            uintPackedOutputStream.write(suffixArray.sequence,1,suffixArray.sequence.length-1);
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
