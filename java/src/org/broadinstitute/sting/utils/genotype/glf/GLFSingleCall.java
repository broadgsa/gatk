package org.broadinstitute.sting.utils.genotype.glf;

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
 * @author aaron
 *         <p/>
 *         Class GLFSingleCall
 *         <p/>
 *         This class represents a single point geneotype call in GLF vernacular
 */
public class GLFSingleCall extends GLFRecord {

    // our likelihoods object
    private double likelihoods[];

    /**
     * create a single
     *
     * @param contig      the contig this record is on
     * @param refBase     the reference base, as a char
     * @param position    the location, as an offset from the start of the contig
     * @param readDepth   the read depth at the specified postion
     * @param rmsMapQ     the root mean square of the mapping quality
     * @param likelihoods the Likelihoods
     */
    public GLFSingleCall(String contig, char refBase, int position, int readDepth, short rmsMapQ, double likelihoods[]) {
        super(contig, refBase, position, (short) GLFRecord.findMin(likelihoods), readDepth, rmsMapQ);
        this.likelihoods = likelihoods;
    }


    /**
     * Write out the record to a binary codec
     *
     * @param out the codec to write to
     */
    void write(BinaryCodec out, GLFRecord lastRec) {
        super.write(out, lastRec);
        short[] adjusted = new short[likelihoods.length];
        // we want to scale our values
        for (int x = 0; x < likelihoods.length; x++) {
            adjusted[x] = GLFRecord.toCappedShort(LIKELIHOOD_SCALE_FACTOR * (Math.round(likelihoods[x]) - this.minimumLikelihood));
        }
        try {
            for (short value : adjusted) {
                out.writeUByte(value);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * return the record type we represent, in this case SINGLE
     *
     * @return RECORD_TYPE.SINGLE
     */
    public RECORD_TYPE getRecordType() {
        return RECORD_TYPE.SINGLE;
    }

    /**
     * return our size in bytes
     *
     * @return number of bytes we represent
     */
    public int getByteSize() {
        return likelihoods.length + super.getByteSize();
    }

    public double[] getLikelihoods() {
        return likelihoods;
    }
}
