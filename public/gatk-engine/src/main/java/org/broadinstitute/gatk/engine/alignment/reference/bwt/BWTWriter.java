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

import org.broadinstitute.gatk.engine.alignment.reference.packing.BasePackedOutputStream;
import org.broadinstitute.gatk.engine.alignment.reference.packing.UnsignedIntPackedOutputStream;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

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
            throw new ReviewedGATKException("Unable to open output file", ex);
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
