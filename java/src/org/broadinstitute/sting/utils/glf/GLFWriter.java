package org.broadinstitute.sting.utils.glf;

import net.sf.samtools.util.BinaryCodec;
import net.sf.samtools.util.BlockCompressedOutputStream;

import java.util.ArrayList;
import java.util.List;
import java.io.File;
import java.io.DataOutputStream;

/**
 *
 * User: aaron
 * Date: May 13, 2009
 * Time: 3:36:18 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 */
public class GLFWriter {
    // our output codec
    private final BinaryCodec outputBinaryCodec;

    public static final short[] glfMagic = {'G', 'L', 'F', '\3'};
    private String headerText = "";
    private String referenceSequenceName = "";
    private long referenceSequenceLength = 0;

    /**
     * The public constructor for creating a GLF object
     *
     * @param headerText            the header text (currently unclear what the contents are)
     * @param referenceSequenceName the reference sequence name
     */
    public GLFWriter( String headerText, String referenceSequenceName, int referenceSequenceLength, File writeTo ) {
        this.headerText = headerText;
        this.referenceSequenceName = referenceSequenceName;
        this.referenceSequenceLength = referenceSequenceLength;
        outputBinaryCodec = new BinaryCodec(new DataOutputStream(new BlockCompressedOutputStream(writeTo)));
        outputBinaryCodec.setOutputFileName(writeTo.toString());
        this.writeHeader();
    }

    /**
     * add a point genotype to the GLF writer
     *
     * @param refBase the reference base, as a char
     * @param genomicLoc the location, as an offset from the previous glf record
     * @param readDepth the read depth at the specified postion
     * @param rmsMapQ the root mean square of the mapping quality
     * @param lhValues the LikelihoodObject, representing the genotype likelyhoods
     */
    public void addPointCall( char refBase,
                              int genomicLoc,
                              int readDepth,
                              short rmsMapQ,
                              LikelihoodObject lhValues ) {

        SinglePointCall call = new SinglePointCall(refBase, genomicLoc,
                readDepth,
                rmsMapQ,
                lhValues.toByteArray(),
                lhValues);
        call.write(this.outputBinaryCodec);
    }

    /**
     * add a variable length (indel, deletion, etc) to the genotype writer
     *
     * @param refBase the reference base
     * @param genomicLoc the location, as an offset from the previous glf record
     * @param readDepth the read depth at the specified postion
     * @param rmsMapQ the root mean square of the mapping quality
     * @param minimumLikelihood the minimum likelihood value
     * @param homozygProb1 the probability of the first homozygous indel allele
     * @param homozygProb2 the probability of the second homozygous indel allele
     * @param heterozygProb the probability of the heterozygote
     * @param indelLength1 the length of the first indel allele
     * @param indelLength2 the length of the second indel allele
     * @param indelSeq1 the sequence for the first indel allele
     * @param indelSeq2 the sequence for the second indel allele
     */
    public void addVarLengthCall( char refBase,
                                  long genomicLoc,
                                  int readDepth,
                                  short rmsMapQ,
                                  short minimumLikelihood,
                                  short homozygProb1,
                                  short homozygProb2,
                                  short heterozygProb,
                                  int indelLength1,
                                  int indelLength2,
                                  char[] indelSeq1,
                                  char[] indelSeq2 ) {

        short[] indexSeqEq1 = new short[indelSeq1.length];
        short[] indexSeqEq2 = new short[indelSeq2.length];
        for (int x = 0; x < indelSeq1.length; x++) {
            indexSeqEq1[x] = new Integer(0x00ff & indelSeq1[x]).shortValue();
        }
        for (int x = 0; x < indelSeq2.length; x++) {
            indexSeqEq2[x] = new Integer(0x00ff & indelSeq2[x]).shortValue();
        }

        VariableLengthCall call = new VariableLengthCall(refBase,
                genomicLoc,
                readDepth,
                minimumLikelihood,
                rmsMapQ,
                homozygProb1,
                homozygProb2,
                heterozygProb,
                indelLength1,
                indelLength2,
                indexSeqEq1,
                indexSeqEq2);

        call.write(this.outputBinaryCodec);

    }

    /**
     * add a GLF record to the output file
     * @param rec the GLF record to write.
     */
    public void addGLFRecord( GLFRecord rec ) {
        rec.write(this.outputBinaryCodec);
    }

    /**
     * Write out the header information for the GLF file
     *
     */
    private void writeHeader() {
        for (int x = 0; x < glfMagic.length; x++) {
            outputBinaryCodec.writeByte(glfMagic[x]);
        }
        if (!(headerText.equals(""))) {
            outputBinaryCodec.writeString(headerText, true, true);
        } else {
            outputBinaryCodec.writeInt(0);
        }
        outputBinaryCodec.writeString(referenceSequenceName, true, true);
        outputBinaryCodec.writeUInt(referenceSequenceLength);
    }

    /**
     * close the file
     *
     */
    public void close() {
        outputBinaryCodec.writeByte((byte) 0);
        outputBinaryCodec.close();
    }

}


