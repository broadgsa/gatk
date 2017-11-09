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

package org.broadinstitute.gatk.tools.walkers.readutils;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class PrintReadsIntegrationTest extends WalkerTest {
    private static class PRTest {
        final String reference;
        final List<String> bam;
        final String args;
        final String md5;

        private PRTest(String reference, String[] bam, String args, String md5) {
            this.reference = reference;
            this.bam = new ArrayList<>();
            this.args = args;
            this.md5 = md5;
            this.bam.addAll(Arrays.asList(bam));
        }

        @Override
        public String toString() {
            return String.format("PRTest(bam='%s', args='%s')", bam, args);
        }
    }

    @DataProvider(name = "PRTest")
    public Object[][] createPrintReadsTestData() {
        return new Object[][]{
                {new PRTest(hg18Reference, new String[]{"HiSeq.1mb.bam"}, "", "83d1454dc01cd2e7458dad4012695f64")},
                {new PRTest(hg18Reference, new String[]{"HiSeq.1mb.bam"}, " -compress 0", "0aec10d19e0dbdfe1d0cbb3eddaf623a")},
                {new PRTest(hg18Reference, new String[]{"HiSeq.1mb.bam"}, " -simplifyBAM", "60255a68df1b8f2fbba373d75274f0de")},
                {new PRTest(hg18Reference, new String[]{"HiSeq.1mb.bam"}, " -n 10", "eb7c6bacca5fee09b8df50880eb81ee6")},
                // See: GATKBAMIndex.getStartOfLastLinearBin(), BAMScheduler.advance(), IntervalOverlapFilteringIterator.advance()
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.bam"}, "", "3d67c398ce2ac1deeddbccbd850380a7")},
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.bam"}, " -L 1", "37cdd8871843693f2650d7b48c8ae1d4")},
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.cram"}, " -L 1:10001 -L GL000192.1:500204", "4d63fd6e977a53e5d9590bd030b40bd0")},
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.bam"}, " -L unmapped", "a834400e3bd69045eb8a9e94131633f5")},
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.bam"}, " -L 1 -L unmapped", "37cdd8871843693f2650d7b48c8ae1d4")},
                {new PRTest(b37KGReference, new String[]{"oneReadAllInsertion.bam"}, "",  "6c04aac25e2136fee395897aac96bea8")},
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam"}, "",  "57a9bc1f7dd4e7717ee796c484bcf45a")},
                // Tests for filtering options
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        "",  "e691d61df10f7614d73c8ecb46c75ee1")},
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        " -readGroup SRR359098",  "7644eab114bf537411218f782d75a6a6")},
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        " -readGroup 20FUK.3 -sn NA12878",  "351d5da29874033e50d29c5c36575a6c")},
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        " -sn na12878",  "9056d852418dd2083f38e3eac1551fcd")},
        };
    }

    @Test(dataProvider = "PRTest")
    public void testPrintReads(PRTest params) {

        StringBuilder inputs = new StringBuilder();
        for (String bam : params.bam) {
            inputs.append(" -I ");
            inputs.append(privateTestDir);
            inputs.append(bam);
        }

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintReads" +
                        " -R " + params.reference +
                        inputs.toString() +
                        params.args +
                        " --no_pg_tag" +
                        " -o %s",
                Collections.singletonList(params.md5));
        executeTest("testPrintReads-"+params.args, spec).getFirst();
    }

    @DataProvider(name = "PRExceptionTest")
    public Object[][] createPrintReadsExceptionTestData() {
        return new Object[][]{
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        "-platform illum",  "")},
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        " -sn NotASample",  "")},
        };
    }

    @Test(dataProvider = "PRExceptionTest")
    public void testPrintReadsException(PRTest params) {

        StringBuilder inputs = new StringBuilder();
        for (String bam : params.bam) {
            inputs.append(" -I ");
            inputs.append(privateTestDir);
            inputs.append(bam);
        }

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintReads" +
                        " -R " + params.reference +
                        inputs.toString() +
                        params.args +
                        " --no_pg_tag" +
                        " -o %s",
                1, UserException.class);
        executeTest("testPrintReadsException-"+params.args, spec);
    }

    @Test
    public void testPrintReadsNoBQSRFile() {

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintReads" +
                        " -R " + hg18Reference +
                        " -I " + privateTestDir + "HiSeq.1mb.bam" +
                        " -BSQR bqsrFile" +
                        " --no_pg_tag" +
                        " -o %s",
                1, UserException.class);
        executeTest("testPrintReadsNoBQSRFile-", spec);
    }

}
