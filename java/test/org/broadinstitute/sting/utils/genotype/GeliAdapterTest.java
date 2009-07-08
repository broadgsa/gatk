package org.broadinstitute.sting.utils.genotype;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.junit.Test;

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
 * 
 * @author aaron 
 * 
 * Class GeliAdapterTest
 *
 * Tests the GeliAdapter class
 */
public class GeliAdapterTest extends BaseTest {


    // private our Geli adapter
    private GenotypeWriter adapter = null;

    /**
     * test out the likelihood object
     */
    @Test
    public void test1() {
        File fl = new File("testFile.txt");
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(2,1,10);
        adapter = new GeliAdapter(fl,header);
        LikelihoodObject obj = new LikelihoodObject(createFakeLikelihoods());
        adapter.addGenotypeCall("chr1",10,100,100,'A',100,obj);
        adapter.close();
    }


    public double[] createFakeLikelihoods() {
        double ret[] = new double[10];
        for (int x = 0; x < 10; x++) {
            ret[x] = (double)(10.0-x) * 10.0;
        }
        return ret;
    }
}
