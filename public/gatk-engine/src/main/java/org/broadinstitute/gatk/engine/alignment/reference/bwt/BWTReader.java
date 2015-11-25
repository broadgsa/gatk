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

import org.broadinstitute.gatk.engine.alignment.reference.packing.BasePackedInputStream;
import org.broadinstitute.gatk.engine.alignment.reference.packing.PackUtils;
import org.broadinstitute.gatk.engine.alignment.reference.packing.UnsignedIntPackedInputStream;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
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
    private FileInputStream inputStream;

    /**
     * Create a new BWT reader.
     * @param inputFile File in which the BWT is stored.
     */
    public BWTReader( File inputFile ) {
        try {
            this.inputStream = new FileInputStream(inputFile);
        }
        catch( FileNotFoundException ex ) {
            throw new ReviewedGATKException("Unable to open input file", ex);
        }
    }

    /**
     * Read a BWT from the input stream.
     * @return The BWT stored in the input stream.
     */
    public BWT read() {
        UnsignedIntPackedInputStream uintPackedInputStream = new UnsignedIntPackedInputStream(inputStream, ByteOrder.LITTLE_ENDIAN);
        BasePackedInputStream basePackedInputStream = new BasePackedInputStream<Integer>(Integer.class, inputStream, ByteOrder.LITTLE_ENDIAN);

        long inverseSA0;
        long[] count;
        SequenceBlock[] sequenceBlocks;

        try {
            inverseSA0 = uintPackedInputStream.read();
            count = new long[PackUtils.ALPHABET_SIZE];
            uintPackedInputStream.read(count);

            long bwtSize = count[PackUtils.ALPHABET_SIZE-1];
            sequenceBlocks = new SequenceBlock[PackUtils.numberOfPartitions(bwtSize,BWT.SEQUENCE_BLOCK_SIZE)];
            
            for( int block = 0; block < sequenceBlocks.length; block++ ) {
                int sequenceStart = block* BWT.SEQUENCE_BLOCK_SIZE;
                int sequenceLength = (int)Math.min(BWT.SEQUENCE_BLOCK_SIZE,bwtSize-sequenceStart);

                long[] occurrences = new long[PackUtils.ALPHABET_SIZE];
                byte[] bwt = new byte[sequenceLength];

                uintPackedInputStream.read(occurrences);
                basePackedInputStream.read(bwt);

                sequenceBlocks[block] = new SequenceBlock(sequenceStart,sequenceLength,new Counts(occurrences,false),bwt);
            }
        }
        catch( IOException ex ) {
            throw new ReviewedGATKException("Unable to read BWT from input stream.", ex);
        }

        return new BWT(inverseSA0, new Counts(count,true), sequenceBlocks);
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
