package org.broadinstitute.sting.utils.glf;

import net.sf.samtools.util.BinaryCodec;


/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * 
 * @author aaron 
 * 
 * Class SinglePointCall
 *
 * This class represents a single point geneotype call in GLF vernacular
 **/
class SinglePointCall extends GLFRecord {

    // our likelyhood array size
    public static final int LIKELYHOOD_SIZE = 10;

    // our record type
    private final RECORD_TYPE type = RECORD_TYPE.SINGLE;

    // our array of likelihoods
    public final short lk[] = new short[LIKELYHOOD_SIZE];

    // our size, we're immutable (the size at least).  In bytes.
    private final int byteSize;


    SinglePointCall(char refBase, int offset, int readDepth, short rmsMapQ, short[] lk, LikelihoodObject minimumLikelihood ) {
        super(refBase,offset,(short)minimumLikelihood.getMinimumValue(),readDepth,rmsMapQ);

        if (lk.length != LIKELYHOOD_SIZE) {
            throw new IllegalArgumentException("SinglePointCall: passed in likelyhood array size != LIKELYHOOD_SIZE");
        }

        System.arraycopy(lk, 0, this.lk, 0, LIKELYHOOD_SIZE);
        byteSize = 9 + lk.length;
    }


    /**
     * Write out the record to a binary codec
     *
     * @param out
     */
    public void write( BinaryCodec out ) {
        super.write(out);
        for (int x = 0; x < LIKELYHOOD_SIZE; x++) {
            out.writeUByte(lk[x]);
        }
    }

    /**
     * return the record type we represent, in this case SINGLE
     * @return RECORD_TYPE.SINGLE
     */
    public RECORD_TYPE getRecordType() {
        return RECORD_TYPE.SINGLE;
    }

    /**
     * return our size in bytes
     * @return number of bytes we represent
     */
    public int getByteSize() {
        return byteSize;
    }

}
