package org.broadinstitute.sting.bwa;

import org.broadinstitute.sting.utils.StingException;

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
    private InputStream inputStream;

    /**
     * Create a new BWT reader.
     * @param inputFile File in which the BWT is stored.
     */
    public BWTReader( File inputFile ) {
        try {
            this.inputStream = new BufferedInputStream(new FileInputStream(inputFile));
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
        int[] occurrences;
        byte[] bwt;

        try {
            inverseSA0 = intPackedInputStream.read();
            occurrences = new int[PackUtils.ALPHABET_SIZE];
            intPackedInputStream.read(occurrences);
            bwt = basePackedInputStream.read(occurrences[PackUtils.ALPHABET_SIZE-1]);
        }
        catch( IOException ex ) {
            throw new StingException("Unable to read BWT from input stream.", ex);
        }

        return new BWT(inverseSA0, occurrences, bwt);
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
