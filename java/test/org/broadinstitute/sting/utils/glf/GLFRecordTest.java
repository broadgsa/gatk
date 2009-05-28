package org.broadinstitute.sting.utils.glf;

import org.junit.Test;
import org.junit.Before;
import net.sf.samtools.util.BinaryCodec;

import java.io.File;


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
public class GLFRecordTest {

    /** some made up values that we use to generate the GLF */
    private final String header = "header";
    private final String referenceSequenceName = "refSeq";
    private final int refLength = 1000;

    private GLFRecord rec;

    @Before
    public void before() {
         rec = new GLFRecord(header, referenceSequenceName, refLength);
        
    }

    /**
     * make a fake snp
     * @param genotype the genotype, 0-15 (AA, AT, AA, ... GG)
     */
    private void addFakeSNP(int genotype, int location) {
        LikelihoodObject obj = new LikelihoodObject();
        obj.setLikelihood(LikelihoodObject.GENOTYPE.values()[genotype],0.5f);
        rec.addSNPCall(location,10,10,obj);
    }


    @Test
    public void basicWrite() {
        File writeTo = new File("testGLF.glf");
        BinaryCodec codec = new BinaryCodec(writeTo, true);
        for (int x = 0; x < 100; x++) {
            addFakeSNP(0,x);
        }
    }

}
