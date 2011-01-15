package org.broadinstitute.sting.utils.genotype.glf;

import net.sf.samtools.util.BinaryCodec;
import org.broadinstitute.sting.utils.Utils;


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
// TODO -- DELETE ME GLF
public class GLFSingleCall extends GLFRecord {

    // our likelihoods array
    private double likelihoods[];
    private double minLikelihood;
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
        super(contig, refBase, position, readDepth, rmsMapQ);
        minLikelihood = GLFRecord.findMin(likelihoods);
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
            adjusted[x] = GLFRecord.toCappedShort(Math.round(LIKELIHOOD_SCALE_FACTOR * (likelihoods[x] - minLikelihood)));
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

    /**
     * this method had to be abstracted so that the underlying records could set the minimum likelihood (ML) in the event
     * that the ML is above 255.  In this case the records need to scale their likelihood values appropriately, and warn the user.
     *
     * @return a short of the minimum likelihood.
     */
    @Override
    protected short calculateMinLikelihood() {
        if (minLikelihood > 255.0) {
            double scale = minLikelihood - 255.0;
            this.minLikelihood = 255.0;
            for (int x = 0; x < this.likelihoods.length; x++)
                this.likelihoods[x] = this.likelihoods[x] - scale;
            Utils.warnUser("GLFRecord: Locus " + this.getContig() + ":" + this.position + " had it's likelihood information scaled, the original likelihood values are unrecoverable");
        }
        return toCappedShort(minLikelihood);
    }

    @Override
    public boolean equals(GLFRecord rec) {
        if (!super.equals(rec)) return false;
        if (!(rec instanceof GLFSingleCall)) return false;
        if (((GLFSingleCall) rec).getLikelihoods().length != this.likelihoods.length) return false;
        for (int x = 0; x < likelihoods.length; x++)
            if (Double.compare(likelihoods[x],((GLFSingleCall) rec).getLikelihoods()[x]) != 0) return false;
        return this.getMinimumLikelihood() == rec.getMinimumLikelihood();
    }

    public double[] getLikelihoods() {
        return likelihoods;
    }


}
