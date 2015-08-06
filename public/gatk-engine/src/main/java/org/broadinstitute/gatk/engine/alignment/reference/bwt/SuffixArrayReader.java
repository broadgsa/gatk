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
import org.broadinstitute.gatk.engine.alignment.reference.packing.UnsignedIntPackedInputStream;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

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
            throw new ReviewedGATKException("Unable to open input file", ex);
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
            throw new ReviewedGATKException("Unable to read BWT from input stream.", ex);
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
            throw new ReviewedGATKException("Unable to close input file", ex);
        }
    }    
}
