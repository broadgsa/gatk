/*
 * Copyright (c) 2011, The Broad Institute
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

        private TestParams(String master, String test, String MD5) {
            super(TestParams.class);
            this.master = new File(master);
            this.test = new File(test);
            this.MD5 = MD5;
        }

        public String toString() {
            return String.format("master=%s,test=%s,md5=%s", master, test, MD5);
        }
    }

    @DataProvider(name = "data")
    public Object[][] createData() {
        new TestParams(testDir + "diffTestMaster.vcf", testDir + "diffTestTest.vcf", "ed377322c615abc7dceb97025076078d");
        new TestParams(testDir + "exampleBAM.bam", testDir + "exampleBAM.simple.bam", "02e46f5d2ebb3d49570850595b3f792e");
        return TestParams.getTests(TestParams.class);
    }

    @Test(enabled = true, dataProvider = "data")
    public void testDiffs(TestParams params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T DiffObjects -R public/testdata/exampleFASTA.fasta "
                        + " -m " + params.master
                        + " -t " + params.test
                        + " -o %s",
                Arrays.asList(params.MD5));
        executeTest("testDiffObjects:"+params, spec).getFirst();
    }
}

