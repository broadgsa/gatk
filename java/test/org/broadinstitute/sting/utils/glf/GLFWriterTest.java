package org.broadinstitute.sting.utils.glf;

import org.junit.Test;
import org.junit.Before;
import org.broadinstitute.sting.BaseTest;
import net.sf.samtools.util.BinaryCodec;
import net.sf.samtools.util.BlockCompressedOutputStream;

import java.io.File;
import java.io.DataOutputStream;
import java.io.IOException;


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
 *         Class GLFRecordTest
 *         <p/>
 *         Tests for the GLFRecord class
 */
public class GLFWriterTest extends BaseTest {

    /** some made up values that we use to generate the GLF */
    private final String header = "";
    private final String referenceSequenceName = "chr1";
    private final int refLength = 1000;
    File writeTo = new File("testGLF.glf");

    private GLFWriter rec;

    @Before
    public void before() {

    }

    /**
     * make a fake snp
     *
     * @param genotype the genotype, 0-15 (AA, AT, AA, ... GG)
     */
    private void addFakeSNP( int genotype, int location ) {
        LikelihoodObject obj = new LikelihoodObject();
        obj.setLikelihood(LikelihoodObject.GENOTYPE.values()[genotype], 128);
        int ran = (int) Math.floor(Math.random() * 4.0);
        char let = 'A';
        switch (ran) {
            case 0:
                let = 'T';
                break;
            case 1:
                let = 'C';
                break;
            case 2:
                let = 'G';
                break;
        }
    try {
        rec.addPointCall(let, location, 10, (short) 10, obj);
    } catch (IllegalArgumentException e) {
        e.printStackTrace();
    }
    }


    @Test
    public void basicWrite() {
        rec = new GLFWriter(header, referenceSequenceName, refLength, writeTo);
        for (int x = 0; x < 100; x++) {
            addFakeSNP((int) Math.round(Math.random() * 9), 1);
        }
        rec.close();
    }

}
