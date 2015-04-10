/*
* Copyright (c) 2012 The Broad Institute
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
                {new PRTest(hg18Reference, new String[]{"HiSeq.1mb.bam"}, "", "5aee1c592f7b0505430df4d4452b8000")},
                {new PRTest(hg18Reference, new String[]{"HiSeq.1mb.bam"}, " -compress 0", "62a542230502c9e54124ebd46242e252")},
                {new PRTest(hg18Reference, new String[]{"HiSeq.1mb.bam"}, " -simplifyBAM", "a054a6618ffa8cd2d1113b005335922b")},
                {new PRTest(hg18Reference, new String[]{"HiSeq.1mb.bam"}, " -n 10", "0e3d1748ad1cb523e3295cab9d09d8fc")},
                // See: GATKBAMIndex.getStartOfLastLinearBin(), BAMScheduler.advance(), IntervalOverlapFilteringIterator.advance()
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.bam"}, "", "d7f23fd77d7dc7cb50d3397f644c6d8a")},
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.bam"}, " -L 1", "c601db95b20248d012b0085347fcb6d1")},
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.bam"}, " -L unmapped", "2d32440e47e8d9d329902fe573ad94ce")},
                {new PRTest(b37KGReference, new String[]{"unmappedFlagReadsInLastLinearBin.bam"}, " -L 1 -L unmapped", "c601db95b20248d012b0085347fcb6d1")},
                {new PRTest(b37KGReference, new String[]{"oneReadAllInsertion.bam"}, "",  "349650b6aa9e574b48a2a62627f37c7d")},
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam"}, "",  "0c1cbe67296637a85e80e7a182f828ab")},
                // Tests for filtering options
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        "",  "b3ae15c8af33fd5badc1a29e089bdaac")},
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        " -readGroup SRR359098",  "8bd867b30539524daa7181efd9835a8f")},
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        " -readGroup 20FUK.3 -sn NA12878",  "93a7bc1b2b1cd27815ed1666cbb4d0cb")},
                {new PRTest(b37KGReference, new String[]{"NA12878.1_10mb_2_10mb.bam", "NA20313.highCoverageRegion.bam"},
                        " -sn na12878",  "52e99cfcf03ff46285d1ba302f8df964")},
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
                Arrays.asList(params.md5));
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
