/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.indels;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.testng.annotations.Test;

import java.util.ArrayList;

public class IndelRealignerLargeScaleTest extends WalkerTest {
    @Test( timeOut = 18000000 )
    public void testHighCoverage() {
        WalkerTestSpec spec = new WalkerTestSpec(

                "-R " + b36KGReference +
                        " -T IndelRealigner" +
                        " -I " + validationDataLocation + "indelRealignerTest.pilot1.veryHighCoverage.bam" +
                        " -L 20:49,500-55,500" +
                        " -o /dev/null" +
                        " -targetIntervals " + validationDataLocation + "indelRealignerTest.pilot1.ceu.intervals",
                 0,
                new ArrayList<String>(0));
        executeTest("testIndelRealignerHighCoverage", spec);
    }

    @Test( timeOut = 18000000 )
    public void testRealigner() {
        WalkerTestSpec spec1 = new WalkerTestSpec(

                "-R " + hg18Reference +
                        " -T IndelRealigner" +
                        " -LOD 5" +
                        " -maxConsensuses 100" +
                        " -greedy 100" +
                        " -known " + GATKDataLocation + "dbsnp_132.hg18.vcf" +
                        " -o /dev/null" +
                        " -I " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.bam" +
                        " -L chr1:1-5,650,000" +
                        " -targetIntervals " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.realigner.intervals",
                 0,
                new ArrayList<String>(0));
        executeTest("testIndelRealignerWholeGenome", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T IndelRealigner" +
                        " -LOD 5" +
                        " -maxConsensuses 100" +
                        " -greedy 100" +
                        " -known " + GATKDataLocation + "dbsnp_132.hg18.vcf" +
                        " -o /dev/null" +
                        " -I " + evaluationDataLocation + "NA12878.ESP.WEx.chr1.bam" +
                        " -L chr1:1-150,000,000" +
                        " -targetIntervals " + evaluationDataLocation + "NA12878.ESP.WEx.chr1.realigner.intervals",
                 0,
                new ArrayList<String>(0));
        executeTest("testIndelRealignerWholeExome", spec2);
    }
}
