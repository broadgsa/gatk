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
                {new PRTest(hg18Reference, new String[]{"HiSeq.1mb.bam"}, "", "0aa3505ba61e05663e629011dd54e423")},
                {new PRTest(hg18Reference, new String[]{"HiSeq.1mb.bam"}, " -compress 0", "0aec10d19e0dbdfe1d0cbb3eddaf623a")},
                {new PRTest(hg18Reference, new String[]{"HiSeq.1mb.bam"}, " -simplifyBAM", "c565d9cd4838a313e7bdb30530c0cf71")},
                {new PRTest(hg18Reference, new String[]{"HiSeq.1mb.bam"}, " -n 10", "917440a38aba707ec0e012168590981a")},
                // See: GATKBAMIndex.getStartOfLastLinearBin(), BAMScheduler.advance(), IntervalOverlapFilteringIterator.advance()
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.bam"}, "", "0b58c903f54e8543a8b2ce1439aa769b")},
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.bam"}, " -L 1", "5b1154cc81dba6bcfe76188e4df8d79c")},
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.cram"}, " -L 1:10001 -L GL000192.1:500204", "e9caf8a0e6ec947cdcbdfc48a4292eb5")},
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.bam"}, " -L unmapped", "cbd3d1d50c8674f79033aa8c36aa3cd1")},
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.bam"}, " -L 1 -L unmapped", "5b1154cc81dba6bcfe76188e4df8d79c")},
                {new PRTest(b37KGReference, new String[]{"oneReadAllInsertion.bam"}, "",  "e212d1799ae797e781b17e630656a9a1")},
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam"}, "",  "0387c61303140d8899fcbfdd3e72ed80")},
                // Tests for filtering options
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        "",  "ad56da66be0bdab5a8992de9617ae6a5")},
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        " -readGroup SRR359098",  "c3bfe28722a665e666098dbb7048a9f1")},
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        " -readGroup 20FUK.3 -sn NA12878",  "8191f8d635d00b1f4d0993b785cc46c5")},
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        " -sn na12878",  "92a85b4223ec45e114f12a1fe6ebbaeb")},
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

}
