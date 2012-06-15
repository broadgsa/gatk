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

package org.broadinstitute.sting.gatk.walkers.diagnostics.targets;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class DiagnoseTargetsIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String singleSample = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
    final String multiSample = validationDataLocation + "CEUTrio.HiSeq.b37.chr20.10_11mb.bam";
    final String L = validationDataLocation + "DT-itest.interval_list";

    private void DTTest(String testName, String args, String md5) {
        String base = String.format("-T DiagnoseTargets -R %s -L %s", REF, L) + " -o %s ";
        WalkerTestSpec spec = new WalkerTestSpec(base + args, Arrays.asList(md5));
        executeTest(testName, spec);
    }

    @Test(enabled = true)
    public void testSingleSample() {
        DTTest("testSingleSample ", "-I " + singleSample + " -max 75", "2df47009571fe83ead779c94be97fe96");
    }

    @Test(enabled = true)
    public void testMultiSample() {
        DTTest("testMultiSample ", "-I " + multiSample, "6f0c070b9671e1d007ce6374c3183014");
    }
}
