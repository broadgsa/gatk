/*
* Copyright 2012-2015 Broad Institute, Inc.
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
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collections;

/**
 * Test the GATK core CRAM parsing mechanism.
 */
public class CramIntegrationTest extends WalkerTest {
    @DataProvider(name="cramData")
    public Object[][] getCRAMData() {
        return new Object[][] {
                {"PrintReads", "exampleBAM.bam", "", "cram", "97470174cd313a4d200b2a96ffd73e99"},
                {"PrintReads", "exampleCRAM.cram", "", "cram", "424c725c4ffe7215e358ecf5abd5e5e8"},
                {"PrintReads", "exampleCRAM.cram", "", "bam", "247805098718dd74b8a871796424d359"},
                {"PrintReads", "exampleCRAM.cram", " -L chr1:200", "bam", "a5b26631cd89f86f6184bcac7bc9c9ca"},
                {"CountLoci", "exampleCRAM.cram", "", "txt", "ade93df31a6150321c1067e749cae9be"},
                {"CountLoci", "exampleCRAM.cram", " -L chr1:200", "txt", "b026324c6904b2a9cb4b88d6d61c81d1"},
                {"CountReads", "exampleCRAM.cram", "", "txt", "4fbafd6948b6529caa2b78e476359875"},
                {"CountReads", "exampleCRAM.cram", " -L chr1:200", "txt", "b026324c6904b2a9cb4b88d6d61c81d1"},
                {"PrintReads", "exampleCRAM.cram", " -L chr1:200 -L chr1:89597", "bam", "24dbd14b60220461f47ec5517962cb7f"},
                {"CountLoci", "exampleCRAM.cram", " -L chr1:200 -L chr1:89597", "txt", "26ab0db90d72e28ad0ba1e22ee510510"},
                {"CountReads", "exampleCRAM.cram", " -L chr1:200 -L chr1:89597", "txt", "6d7fce9fee471194aa8b5b6e47267f03"},
                {"PrintReads", "exampleCRAM-nobai-withcrai.cram", " -L chr1:200 -L chr1:89597", "bam", "84bee5063d8fa0d07e7c3ff7e825ae3a"},
                {"CountLoci", "exampleCRAM-nobai-withcrai.cram", " -L chr1:200 -L chr1:89597", "txt", "26ab0db90d72e28ad0ba1e22ee510510"},
                {"CountReads", "exampleCRAM-nobai-withcrai.cram", " -L chr1:200 -L chr1:89597", "txt", "6d7fce9fee471194aa8b5b6e47267f03"},
        };
    }

    @Test(dataProvider = "cramData")
    public void testCram(String walker, String input, String args, String ext, String md5) {
        WalkerTestSpec spec = new WalkerTestSpec(
                " -T Test" + walker + "Walker" +
                    " -I " + publicTestDir + input +
                    " -R " + exampleFASTA +
                    args +
                    " -o %s",
                1, // just one output file
                Collections.singletonList(ext),
                Collections.singletonList(md5));
        executeTest(String.format("testCram %s %s -> %s: %s", walker, input, ext, args), spec);
    }

    @DataProvider(name = "cramNoIndexData")
    public Object[][] getCramNoIndexData() {
        return new Object[][]{
                {"exampleCRAM-nobai-nocrai.cram"},
        };
    }

    @Test(dataProvider = "cramNoIndexData")
    public void testCramNoIndex(String input) {
        WalkerTestSpec spec = new WalkerTestSpec(
                " -T TestPrintReadsWalker" +
                        " -I " + publicTestDir + input +
                        " -R " + exampleFASTA,
                0,
                UserException.class);
        executeTest(String.format("testCramNoIndex %s", input), spec);
    }
}
