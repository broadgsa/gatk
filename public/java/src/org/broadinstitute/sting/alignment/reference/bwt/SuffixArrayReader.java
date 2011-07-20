package org.broadinstitute.sting.alignment.reference.bwt;

import org.broadinstitute.sting.alignment.reference.packing.PackUtils;
import org.broadinstitute.sting.alignment.reference.packing.UnsignedIntPackedInputStream;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
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
    private FileInputStream inputStream;

    /**
     * BWT to use to fill in missing data.
     */
    private BWT bwt;

    /**
     * Create a new suffix array reader.
     * @param inputFile File in which the suffix array is stored.
     * @param bwt BWT to use when filling in missing data.
     */
    public SuffixArrayReader(File inputFile, BWT bwt) {
        try {
            this.inputStream = new FileInputStream(inputFile);
            this.bwt = bwt;
        }
        catch( FileNotFoundException ex ) {
            throw new ReviewedStingException("Unable to open input file", ex);
        }
    }

    /**
     * Read a suffix array from the input stream.
     * @return The suffix array stored in the input stream.
     */
    public SuffixArray read() {
        UnsignedIntPackedInputStream uintPackedInputStream = new UnsignedIntPackedInputStream(inputStream, ByteOrder.LITTLE_ENDIAN);

        long inverseSA0;
        long[] occurrences;
        long[] suffixArray;
        int suffixArrayInterval;

        try {
            inverseSA0 = uintPackedInputStream.read();
            occurrences = new long[PackUtils.ALPHABET_SIZE];
            uintPackedInputStream.read(occurrences);
            // Throw away the suffix array size in bytes and use the occurrences table directly.
            suffixArrayInterval = (int)uintPackedInputStream.read();
            suffixArray = new long[(int)((occurrences[occurrences.length-1]+suffixArrayInterval-1)/suffixArrayInterval)];
            uintPackedInputStream.read(suffixArray);
        }
        catch( IOException ex ) {
            throw new ReviewedStingException("Unable to read BWT from input stream.", ex);
        }

        return new SuffixArray(inverseSA0, new Counts(occurrences,true), suffixArray, suffixArrayInterval, bwt);
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
