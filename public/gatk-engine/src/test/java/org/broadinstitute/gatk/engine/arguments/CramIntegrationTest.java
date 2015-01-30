/*
* Copyright (c) 2012 The Broad Institute
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.arguments;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Test the GATK core CRAM parsing mechanism.
 */
public class CramIntegrationTest extends WalkerTest {
    @DataProvider(name="cramData")
    public Object[][] getCRAMData() {
        return new Object[][] {
                {"PrintReads", "exampleBAM.bam", "", "cram", "026ebc00c2a8f9832e37f1a6a0f53521"},
                //{"PrintReads", "exampleCRAM.cram", "", "cram", "026ebc00c2a8f9832e37f1a6a0f53521"}, https://github.com/samtools/htsjdk/issues/148
                {"PrintReads", "exampleCRAM.cram", "", "bam", "99e5f740b43594a5b8e5bc1a007719e0"},
                {"PrintReads", "exampleCRAM-noindex.cram", "", "bam", "99e5f740b43594a5b8e5bc1a007719e0"},
                {"PrintReads", "exampleCRAM.cram", " -L chr1:200", "bam", "072435e8272411c31b2234f851706384"},
                {"PrintReads", "exampleCRAM-noindex.cram", " -L chr1:200", "bam", "072435e8272411c31b2234f851706384"},
                {"CountLoci", "exampleCRAM.cram", "", "txt", "ade93df31a6150321c1067e749cae9be"},
                {"CountLoci", "exampleCRAM-noindex.cram", "", "txt", "ade93df31a6150321c1067e749cae9be"},
                {"CountLoci", "exampleCRAM.cram", " -L chr1:200", "txt", "b026324c6904b2a9cb4b88d6d61c81d1"},
                {"CountLoci", "exampleCRAM-noindex.cram", " -L chr1:200", "txt", "b026324c6904b2a9cb4b88d6d61c81d1"},
                {"CountReads", "exampleCRAM.cram", "", "txt", "4fbafd6948b6529caa2b78e476359875"},
                {"CountReads", "exampleCRAM-noindex.cram", "", "txt", "4fbafd6948b6529caa2b78e476359875"},
                {"CountReads", "exampleCRAM.cram", " -L chr1:200", "txt", "b026324c6904b2a9cb4b88d6d61c81d1"},
                {"CountReads", "exampleCRAM-noindex.cram", " -L chr1:200", "txt", "b026324c6904b2a9cb4b88d6d61c81d1"},
                {"PrintReads", "exampleCRAM.cram", " -L chr1:200 -L chr1:89597", "bam", "9598062587ad8d2ec596a8ecb19be979"},
                {"CountLoci", "exampleCRAM.cram", " -L chr1:200 -L chr1:89597", "txt", "26ab0db90d72e28ad0ba1e22ee510510"},
                {"CountReads", "exampleCRAM.cram", " -L chr1:200 -L chr1:89597", "txt", "6d7fce9fee471194aa8b5b6e47267f03"},
        };
    }

    @Test(dataProvider = "cramData")
    public void testCRAM(String walker, String input, String args, String ext, String md5) {
        WalkerTestSpec spec = new WalkerTestSpec(
                " -T Test" + walker + "Walker" +
                    " -I " + publicTestDir + input +
                    " -R " + exampleFASTA +
                    args +
                    " -o %s",
                1, // just one output file
                Arrays.asList(ext),
                Arrays.asList(md5));
        executeTest(String.format("testCRAM %s %s -> %s: %s", walker, input, ext, args), spec);
    }
}
