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

package org.broadinstitute.gatk.tools.walkers.readutils;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.testng.annotations.Test;

import java.util.ArrayList;

public class PrintReadsLargeScaleTest extends WalkerTest {
    @Test( timeOut = 1000 * 60 * 60 * 20 ) // 60 seconds * 60 seconds / minute * 20 minutes
    public void testRealignerTargetCreator() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + b37KGReference +
                        " -T PrintReads" +
                        " -I " + evaluationDataLocation + "CEUTrio.HiSeq.WEx.b37.NA12892.clean.dedup.recal.1.bam" +
                        " -o /dev/null",
                 0,
                new ArrayList<String>(0));
        executeTest("testPrintReadsWholeExomeChr1", spec);
    }
}
