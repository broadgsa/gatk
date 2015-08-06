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

import org.broadinstitute.gatk.engine.alignment.reference.packing.UnsignedIntPackedOutputStream;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

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
            throw new ReviewedGATKException("Unable to open input file", ex);
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
            throw new ReviewedGATKException("Unable to read BWT from input stream.", ex);
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
            throw new ReviewedGATKException("Unable to close input file", ex);
        }
    }
}
