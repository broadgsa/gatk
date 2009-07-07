package org.broadinstitute.sting.utils.genotype.glf;

import net.sf.samtools.util.BinaryCodec;
import net.sf.samtools.util.BlockCompressedOutputStream;

import java.io.File;
import java.io.DataOutputStream;
import java.util.Iterator;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.LikelihoodObject;

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
 * an object for reading in GLF files
 */
public class GLFReader implements Iterator<GLFRecord> {

    // our next record
    private GLFRecord nextRecord = null;

    // the glf magic number, which identifies a properly formatted GLF file
    public static final short[] glfMagic = {'G', 'L', 'F', '\3'};

    // our input codec
    private final BinaryCodec inputBinaryCodec;

    // our header string
    private String headerStr;

    // our reference name
    private String referenceName;

    // reference length
    private int referenceLength;

    GLFReader( File readFrom ) {
        inputBinaryCodec = new BinaryCodec(new DataOutputStream(new BlockCompressedOutputStream(readFrom)));
        inputBinaryCodec.setInputFileName(readFrom.getName());

        // first verify that it's a valid GLF
        for (short s: glfMagic) {
            if (inputBinaryCodec.readUByte() != s) throw new StingException("Verification of GLF format failed: magic string doesn't match)");
        }

        // get the header string
        headerStr = inputBinaryCodec.readLengthAndString(false);

        // get the reference name
        referenceName = inputBinaryCodec.readLengthAndString(true);

        // get the reference length - this may be a problem storing an unsigned int into a signed int.  but screw it.
        referenceLength = (int)inputBinaryCodec.readUInt();

        // get the next record
        nextRecord = next();
    }

    private SinglePointCall generateSPC(char refBase, BinaryCodec inputBinaryCodec) {
        int offset = (int)inputBinaryCodec.readUInt();
        long depth = inputBinaryCodec.readUInt();
        short min_lk = (short)((depth & 0x00000000ff000000) >> 24);
        int readDepth = (int)(depth & 0x0000000000ffffff);
        short rmsMapping = inputBinaryCodec.readUByte();
        double[] lkValues = new double[LikelihoodObject.GENOTYPE.values().length];
        for (int x = 0; x < LikelihoodObject.GENOTYPE.values().length; x++) {
            lkValues[x] = inputBinaryCodec.readUByte();
        }
        return new SinglePointCall(refBase,offset,readDepth,rmsMapping,lkValues);
    }


    private VariableLengthCall generateVLC(char refBase, BinaryCodec inputBinaryCodec) {
        int offset = (int)inputBinaryCodec.readUInt();
        int depth = (int)inputBinaryCodec.readUInt();
        short min_lk = (short)((depth & 0x00000000ff000000) >> 24);
        int readDepth = (depth & 0x0000000000ffffff);
        short rmsMapping = inputBinaryCodec.readUByte();
        short lkHom1 = inputBinaryCodec.readUByte();
        short lkHom2 = inputBinaryCodec.readUByte();
        short lkHet = inputBinaryCodec.readUByte();
        int indelLen1 = (int)inputBinaryCodec.readShort();
        int indelLen2 = (int)inputBinaryCodec.readShort();
        short[] indelSeq1 = new short[indelLen1];
        short[] indelSeq2 = new short[indelLen2];
        for (int x = 0; x < indelLen1; x++) {
            indelSeq1[x] = inputBinaryCodec.readUByte();
        }
        for (int x = 0; x < indelLen2; x++) {
            indelSeq2[x] = inputBinaryCodec.readUByte();    
        }
        return new VariableLengthCall(refBase,offset,readDepth,rmsMapping,lkHom1,lkHom2,lkHet,indelSeq1,indelSeq2);
    }

    @Override
    public boolean hasNext() {
        return (nextRecord != null);
    }

    @Override
    public GLFRecord next() {
        short firstBase = inputBinaryCodec.readUByte();
        byte recordType = (byte)(firstBase & 0x00f0 >> 4);
        char refBase = (char)(firstBase & 0x000f);

        GLFRecord ret = nextRecord;
        if (recordType == 1) {
            nextRecord = generateSPC(refBase, inputBinaryCodec);
        }
        else if (recordType == 2) {
            nextRecord =  generateVLC(refBase, inputBinaryCodec);
        }
        else if (recordType == 0){
            nextRecord = null;
        }
        return ret;
    }

    @Override
    public void remove() {
        throw new StingException("I'm Sorry Dave, I can't let you do that (also GLFReader doesn't support remove()).");
    }
}
