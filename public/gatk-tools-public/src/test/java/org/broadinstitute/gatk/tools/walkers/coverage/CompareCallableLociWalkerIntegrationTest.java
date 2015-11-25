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

package org.broadinstitute.gatk.tools.walkers.coverage;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class CompareCallableLociWalkerIntegrationTest extends WalkerTest {
    final static String commonArgs = "-R " + hg18Reference + " -T CompareCallableLoci --comp1:Bed " + validationDataLocation + "1kg_slx.chr1_10mb.callable.bed --comp2:Bed " + validationDataLocation + "ga2_slx.chr1_10mb.callable.bed -o %s";

    @Test
    public void testCompareCallableLociWalker1() {
        String gatk_args = commonArgs + " -L chr1:1-10,000,000";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 1, Arrays.asList("70efebac55ff210bb022e9e22ed80a95"));
        executeTest("CompareCallableLoci Walker", spec);
    }
}
