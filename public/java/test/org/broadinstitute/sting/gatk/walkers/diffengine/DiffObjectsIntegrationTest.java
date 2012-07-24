/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.diffengine;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class DiffObjectsIntegrationTest extends WalkerTest {
    private class TestParams extends TestDataProvider {
        public File master, test;
        public String MD5;
        public boolean doPairwise;

        private TestParams(String master, String test, final boolean doPairwise, String MD5) {
            super(TestParams.class);
            this.master = new File(master);
            this.test = new File(test);
            this.MD5 = MD5;
            this.doPairwise = doPairwise;
        }

        public String toString() {
            return String.format("master=%s,test=%s,md5=%s", master, test, MD5);
        }
    }

    @DataProvider(name = "data")
    public Object[][] createData() {
        new TestParams(privateTestDir + "diffTestMaster.vcf", privateTestDir + "diffTestTest.vcf", true, "aea3d5df32a2acd400da48d06b4dbc60");
        new TestParams(publicTestDir + "exampleBAM.bam", publicTestDir + "exampleBAM.simple.bam", true, "3f46f5a964f7c34015d972256fe49a35");
        new TestParams(privateTestDir + "diffTestMaster.vcf", privateTestDir + "diffTestTest.vcf", false, "e71e23e7ebfbe768e59527bc62f8918d");
        new TestParams(publicTestDir + "exampleBAM.bam", publicTestDir + "exampleBAM.simple.bam", false, "47bf16c27c9e2c657a7e1d13f20880c9");
        return TestParams.getTests(TestParams.class);
    }

    @Test(enabled = true, dataProvider = "data")
    public void testDiffs(TestParams params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T DiffObjects -R public/testdata/exampleFASTA.fasta "
                        + " -m " + params.master
                        + " -t " + params.test
                        + (params.doPairwise ? " -doPairwise " : "")
                        + " -o %s",
                Arrays.asList(params.MD5));
        executeTest("testDiffObjects:"+params, spec).getFirst();
    }
}

