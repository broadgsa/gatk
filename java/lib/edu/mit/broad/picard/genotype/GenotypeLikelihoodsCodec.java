/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.genotype;

import java.io.InputStream;
import java.io.OutputStream;

import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.sam.util.RuntimeEOFException;
import edu.mit.broad.sam.util.SortingCollection;

public class GenotypeLikelihoodsCodec implements SortingCollection.Codec<GenotypeLikelihoods> {
    private static final int SIG_FIG_MULTIPLIER = 100;
    private static final short BLOCK_SIZE = 12 + 10 * 4;

    private OutputStream os;
    private InputStream is;
    private BinaryCodec binaryCodec;

    /** Returns a new genotype likelihood codec. */
    public SortingCollection.Codec<GenotypeLikelihoods> clone() {
        return new GenotypeLikelihoodsCodec();
    }

    /**
     * Write object to OutputStream.
     *
     * @param genotypeLikelihoods what to write
     */
    public void encode(final GenotypeLikelihoods genotypeLikelihoods) {
        this.binaryCodec.writeShort(BLOCK_SIZE);
        this.binaryCodec.writeUInt(genotypeLikelihoods.getReferenceIndex());
        this.binaryCodec.writeUInt(genotypeLikelihoods.getPosition());
        this.binaryCodec.writeByte(genotypeLikelihoods.getReferenceBase());
        this.binaryCodec.writeUShort(genotypeLikelihoods.getNumReads());
        this.binaryCodec.writeByte(genotypeLikelihoods.getMaxMappingQuality());
        
        for (int i = 0; i < genotypeLikelihoods.getLikelihoods().length; i++) {
            writeLikelihood(genotypeLikelihoods.getLikelihoods()[i]);
        }
    }

    /**
     * Read the next record from the input stream and convert into a java object.
     *
     * @return null if no more records.  Should throw exception if EOF is encountered in the middle of
     *         a record.
     */
    public GenotypeLikelihoods decode() {
        int recordLength = 0;
        try {
            recordLength = this.binaryCodec.readShort();
        } catch (RuntimeEOFException e) {
            return null;
        }
        if (recordLength != BLOCK_SIZE) {
            throw new GeliException("Invalid record length: " + recordLength);
        }
        
        final GenotypeLikelihoods genotypeLikelihoods = new GenotypeLikelihoods();
        genotypeLikelihoods.setReferenceIndex(this.binaryCodec.readUInt());
        genotypeLikelihoods.setPosition(this.binaryCodec.readUInt());
        genotypeLikelihoods.setReferenceBase(this.binaryCodec.readByte());
        genotypeLikelihoods.setNumReads(this.binaryCodec.readUShort());
        genotypeLikelihoods.setMaxMappingQuality(this.binaryCodec.readByte());
        
        int bestIndex = -1;
        int secondBestIndex = -1;
        for (int i = 0; i < genotypeLikelihoods.getLikelihoods().length; i++) {
            float likelihood = readLikelihood();
            genotypeLikelihoods.getLikelihoods()[i] = likelihood;
            
            if (bestIndex == -1 || genotypeLikelihoods.getLikelihood(bestIndex) < likelihood) {
                secondBestIndex = bestIndex;
                bestIndex = i;
            } else if (secondBestIndex == -1 || genotypeLikelihoods.getLikelihood(secondBestIndex) < likelihood) {
                secondBestIndex = i;
            }
        }
        genotypeLikelihoods.setBestLikelihoodIndex(bestIndex);
        genotypeLikelihoods.setSecondBestLikelihoodIndex(secondBestIndex);
        
        return genotypeLikelihoods;
    }

    /**
     * Where to write encoded output
     *
     * @param os
     */
    public void setOutputStream(final OutputStream os) {
        this.os = os;
        this.binaryCodec = new BinaryCodec(os);
    }

    /**
     * Where to read encoded input from
     *
     * @param is
     */
    public void setInputStream(final InputStream is) {
        this.is = is;
        this.binaryCodec = new BinaryCodec(is);
    }
    
    private void writeLikelihood(float likelihood) {
        float shiftedLikelihood = likelihood * SIG_FIG_MULTIPLIER;
        this.binaryCodec.writeInt((int) Math.round(shiftedLikelihood));
    }

    /**
     * @return
     */
    private float readLikelihood() {
        float likelihood = (float) this.binaryCodec.readInt() / SIG_FIG_MULTIPLIER;
        return likelihood;
    }

}
