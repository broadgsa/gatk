package org.broadinstitute.sting.bwa;

import org.broadinstitute.sting.utils.StingException;

import java.io.*;
import java.nio.ByteOrder;

/**
 * A reader for suffix arrays in permanent storage.
 *
 * @author mhanna
 * @version 0.1
 */
public class SuffixArrayReader {
    /**
     * Input stream from which to read suffix array data.
     */
    private InputStream inputStream;

    /**
     * Create a new suffix array reader.
     * @param inputFile File in which the suffix array is stored.
     */
    public SuffixArrayReader( File inputFile ) {
        try {
            this.inputStream = new BufferedInputStream(new FileInputStream(inputFile));
        }
        catch( FileNotFoundException ex ) {
            throw new StingException("Unable to open input file", ex);
        }
    }

    /**
     * Read a suffix array from the input stream.
     * @return The suffix array stored in the input stream.
     */
    public SuffixArray read() {
        IntPackedInputStream intPackedInputStream = new IntPackedInputStream(inputStream, ByteOrder.LITTLE_ENDIAN);

        int inverseSA0;
        int[] occurrences;
        int[] suffixArray;

        try {
            inverseSA0 = intPackedInputStream.read();
            occurrences = new int[PackUtils.ALPHABET_SIZE];
            intPackedInputStream.read(occurrences);
            // Throw away the suffix array size in bytes and use the occurrences table directly.
            intPackedInputStream.read();
            int suffixArraySize = occurrences[PackUtils.ALPHABET_SIZE-1]+1;
            suffixArray = new int[suffixArraySize];
            intPackedInputStream.read(suffixArray);
        }
        catch( IOException ex ) {
            throw new StingException("Unable to read BWT from input stream.", ex);
        }

        return new SuffixArray(inverseSA0, new Counts(occurrences,true), suffixArray);
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
