package org.broadinstitute.sting.utils.genotype.glf;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.BinaryCodec;
import net.sf.samtools.util.BlockCompressedOutputStream;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.genotype.*;

import java.io.DataOutputStream;
import java.io.File;
import java.io.OutputStream;
import java.util.List;
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
 * @version 1.0
 *          <p/>
 *          This class writes GLF files. You can either specify GLFRecords, or programaticly generate
 *          single and variable length genotype calls using the provided functions.  When you've finished
 *          generating GLF records, make sure you close the file.
 */
public class GLFWriter implements GenotypeWriter {
    // our output codec
    private final BinaryCodec outputBinaryCodec;

    // the glf magic number, which identifies a properly formatted GLF file
    public static final short[] glfMagic = {'G', 'L', 'F', '\3'};

    // our header text, reference sequence name (i.e. chr1), and it's length
    private String headerText = "";
    private String referenceSequenceName = null;
    private long referenceSequenceLength = 0;

    // the last position written
    private int lastPos = 1;

    /**
     * The public constructor for creating a GLF object
     *
     * @param headerText the header text (currently unclear what the contents are)
     * @param writeTo    the location to write to
     */
    public GLFWriter(String headerText, File writeTo) {
        this.headerText = headerText;
        outputBinaryCodec = new BinaryCodec(new DataOutputStream(new BlockCompressedOutputStream(writeTo)));
        outputBinaryCodec.setOutputFileName(writeTo.toString());
        this.writeHeader();
    }

    /**
     * The public constructor for creating a GLF object
     *
     * @param headerText the header text (currently unclear what the contents are)
     * @param writeTo    the location to write to
     */
    public GLFWriter(String headerText, OutputStream writeTo) {
        this.headerText = headerText;
        outputBinaryCodec = new BinaryCodec(writeTo);
        outputBinaryCodec.setOutputFileName(writeTo.toString());
        this.writeHeader();
    }


    /**
     * add a point genotype to the GLF writer
     *
     * @param contig     the name of the contig you're calling in
     * @param refBase    the reference base, as a char
     * @param genomicLoc the location the location on the reference contig
     * @param readDepth  the read depth at the specified postion
     * @param rmsMapQ    the root mean square of the mapping quality
     * @param lhValues   the GenotypeLikelihoods object, representing the genotype likelyhoods
     */
    public void addGenotypeCall(SAMSequenceRecord contig,
                                int genomicLoc,
                                float rmsMapQ,
                                char refBase,
                                int readDepth,
                                LikelihoodObject lhValues) {

        // check if we've jumped to a new contig
        checkSequence(contig.getSequenceName(), contig.getSequenceLength());

        SinglePointCall call = new SinglePointCall(refBase,
                                                   genomicLoc - lastPos,
                                                   readDepth,
                                                   (short) rmsMapQ,
                                                   lhValues.toDoubleArray());
        lastPos = genomicLoc;
        call.write(this.outputBinaryCodec);
    }

    /**
     * Add a genotype, given a genotype call
     *
     * @param call the genotype call
     */
    public void addGenotypeCall(Genotype call) {
        if ( !(call instanceof GLFGenotypeCall) )
            throw new IllegalArgumentException("Only GeliGenotypeCalls should be passed in to the Geli writers");
        GLFGenotypeCall gCall = (GLFGenotypeCall)call;

        char ref = gCall.getReference();

        // get likelihood information if available
        LikelihoodObject obj = new LikelihoodObject(gCall.getLikelihoods(), LikelihoodObject.LIKELIHOOD_TYPE.LOG);
        obj.setLikelihoodType(LikelihoodObject.LIKELIHOOD_TYPE.NEGATIVE_LOG);  // transform! ... to negitive log likelihoods

        // calculate the RMS mapping qualities and the read depth
        double rms = 0.0;
        if ( gCall.getPileup() != null )
            rms = calculateRMS(gCall.getPileup().getReads());
        int readCount = gCall.getReadCount();
        this.addGenotypeCall(GenomeLocParser.getContigInfo(gCall.getLocation().getContig()),(int)gCall.getLocation().getStart(),(float)rms,ref,readCount,obj);
    }


    /**
     * calculate the rms , given the read pileup
     *
     * @param reads the read array
     *
     * @return the rms of the read mapping qualities
     */
    private double calculateRMS(List<SAMRecord> reads) {
        int[] qualities = new int[reads.size()];
        for (int i = 0; i < reads.size(); i++) {
            qualities[i] = reads.get(i).getMappingQuality();
        }
        return MathUtils.rms(qualities);
    }
    /**
     * add a variable length (indel, deletion, etc) to the genotype writer
     *
     * @param contig        the name of the contig you're calling in
     * @param refBase       the reference base
     * @param genomicLoc    the location on the reference contig
     * @param readDepth     the read depth at the specified postion
     * @param rmsMapQ       the root mean square of the mapping quality
     * @param firstHomZyg   the first homozygous call
     * @param secondHomZyg  the second homozygous call
     * @param hetLikelihood the negitive log likelihood of the heterozygote,  from 0 to 255
     */
    public void addVariableLengthCall(SAMSequenceRecord contig,
                                      int genomicLoc,
                                      float rmsMapQ,
                                      int readDepth,
                                      char refBase,
                                      IndelLikelihood firstHomZyg,
                                      IndelLikelihood secondHomZyg,
                                      byte hetLikelihood) {

        // check if we've jumped to a new contig
        checkSequence(contig.getSequenceName(), contig.getSequenceLength());

        // normalize the two
        VariableLengthCall call = new VariableLengthCall(refBase,
                                                         genomicLoc - lastPos,
                                                         readDepth,
                                                         (short) rmsMapQ,
                                                         firstHomZyg.getLikelihood(),
                                                         secondHomZyg.getLikelihood(),
                                                         hetLikelihood,
                                                         firstHomZyg.getLengthOfIndel(),
                                                         firstHomZyg.getIndelSequence(),
                                                         secondHomZyg.getLengthOfIndel(),
                                                         secondHomZyg.getIndelSequence());
        lastPos = genomicLoc;
        call.write(this.outputBinaryCodec);
    }

    /**
     * add a no call to the genotype file, if supported.
     *
     * @param position the position
     */
    public void addNoCall(int position) {
        // glf doesn't support this operation
        throw new UnsupportedOperationException("GLF doesn't support a 'no call' call.");
    }

    /**
     * add a GLF record to the output file
     *
     * @param contigName   the contig name
     * @param contigLength the contig length
     * @param rec          the GLF record to write.
     */
    public void addGLFRecord(String contigName, int contigLength, GLFRecord rec) {
        checkSequence(contigName, contigLength);
        rec.write(this.outputBinaryCodec);
    }

    /**
     * Write out the header information for the GLF file.  The header contains
     * the magic number, the length of the header text, the text itself, the reference
     * sequence (null terminated) preceeded by it's length, and the the genomic
     * length of the reference sequence.
     */
    private void writeHeader() {
        for (int x = 0; x < glfMagic.length; x++) {
            outputBinaryCodec.writeUByte(glfMagic[x]);
        }
        if (!(headerText.equals(""))) {
            outputBinaryCodec.writeString(headerText, true, true);
        } else {
            outputBinaryCodec.writeInt(0);
        }
    }

    /**
     * check to see if we've jumped to a new contig
     *
     * @param sequenceName the name for the sequence
     * @param seqLength    the sequence length
     */
    private void checkSequence(String sequenceName, int seqLength) {
        if ((referenceSequenceName == null) || (!referenceSequenceName.equals(sequenceName))) {
            if (this.referenceSequenceName != null) { // don't write the record the first time
                this.writeEndRecord();
            }
            referenceSequenceName = sequenceName;
            referenceSequenceLength = seqLength;
            lastPos = 1;
            addSequence();
        }
    }


    /** add a sequence definition to the glf */
    private void addSequence() {
        outputBinaryCodec.writeString(referenceSequenceName, true, true);
        outputBinaryCodec.writeUInt(referenceSequenceLength);
    }

    /** write end record */
    private void writeEndRecord() {
        outputBinaryCodec.writeUByte((short) 0);
    }


    /**
     * close the file.  You must close the file to ensure any remaining data gets
     * written out.
     */
    public void close() {
        writeEndRecord();
        outputBinaryCodec.close();
    }

    /**
     * get the reference sequence
     *
     * @return the reference sequence
     */
    public String getReferenceSequenceName() {
        return referenceSequenceName;
    }


    /**
     * add a multi-sample call if we support it
     *
     * @param genotypes the list of genotypes
     */
    public void addMultiSampleCall(List<Genotype> genotypes, VariationCall metadata) {
        throw new UnsupportedOperationException("GLF writer doesn't support multisample calls");
    }

    /** @return true if we support multisample, false otherwise */
    public boolean supportsMultiSample() {
        return false;
    }
}


