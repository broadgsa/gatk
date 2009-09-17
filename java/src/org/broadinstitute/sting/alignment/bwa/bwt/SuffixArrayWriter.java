package org.broadinstitute.sting.alignment.bwa.bwt;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.alignment.bwa.packing.IntPackedOutputStream;

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
     * @param inputFile File in which the suffix array is stored.
     */
    public SuffixArrayWriter( File inputFile ) {
        try {
            this.outputStream = new BufferedOutputStream(new FileOutputStream(inputFile));
        }
        catch( FileNotFoundException ex ) {
            throw new StingException("Unable to open input file", ex);
        }
    }

    /**
     * Write a suffix array to the output stream.
     * @param suffixArray suffix array to write.
     */
    public void write(SuffixArray suffixArray) {
        IntPackedOutputStream intPackedOutputStream = new IntPackedOutputStream(outputStream, ByteOrder.LITTLE_ENDIAN);

        try {
            intPackedOutputStream.write(suffixArray.inverseSA0);
            intPackedOutputStream.write(suffixArray.occurrences.toArray(true));
            // How frequently the suffix array entry is placed.
            intPackedOutputStream.write(1);
            // Length of the suffix array.
            intPackedOutputStream.write(suffixArray.length()-1);
            intPackedOutputStream.write(suffixArray.sequence, 1, suffixArray.length()-1);
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
