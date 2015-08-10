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

package org.broadinstitute.gatk.tools.walkers.diffengine;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.BaseTest;
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
        new TestParams(privateTestDir + "diffTestMaster.vcf", privateTestDir + "diffTestTest.vcf", true, "71869ddf9665773a842a9def4cc5f3c8");
        new TestParams(publicTestDir + "exampleBAM.bam", publicTestDir + "exampleBAM.simple.bam", true, "cec7c644c84ef9c96aacaed604d9ec9b");
        new TestParams(privateTestDir + "diffTestMaster.vcf", privateTestDir + "diffTestTest.vcf", false, "47546e03344103020e49d8037a7e0727");
        new TestParams(publicTestDir + "exampleBAM.bam", publicTestDir + "exampleBAM.simple.bam", false, "d27b37f7a366c8dacca5cd2590d3c6ce");
        return TestParams.getTests(TestParams.class);
    }

    @Test(enabled = true, dataProvider = "data")
    public void testDiffs(TestParams params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T DiffObjects -R " + publicTestDir + "exampleFASTA.fasta "
                        + " -m " + params.master
                        + " -t " + params.test
                        + (params.doPairwise ? " -doPairwise " : "")
                        + " -o %s",
                Arrays.asList(params.MD5));
        executeTest("testDiffObjects:"+params, spec).getFirst();
    }
}

