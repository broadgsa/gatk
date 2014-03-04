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

package org.broadinstitute.sting.gatk.walkers.readutils;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public class PrintReadsIntegrationTest extends WalkerTest {
    private static class PRTest {
        final String reference;
        final String bam;
        final String args;
        final String md5;

        private PRTest(String reference, String bam, String args, String md5) {
            this.reference = reference;
            this.bam = bam;
            this.args = args;
            this.md5 = md5;
        }

        @Override
        public String toString() {
            return String.format("PRTest(bam='%s', args='%s')", bam, args);
        }
    }

    @DataProvider(name = "PRTest")
    public Object[][] createPrintReadsTestData() {
        return new Object[][]{
                {new PRTest(hg18Reference, "HiSeq.1mb.bam", "", "fa9c66f66299fe5405512ac36ec9d0f2")},
                {new PRTest(hg18Reference, "HiSeq.1mb.bam", " -compress 0", "488eb22abc31c6af7cbb1a3d41da1507")},
                {new PRTest(hg18Reference, "HiSeq.1mb.bam", " -simplifyBAM", "1510dc4429f3ed49caf96da41e8ed396")},
                {new PRTest(hg18Reference, "HiSeq.1mb.bam", " -n 10", "0e3d1748ad1cb523e3295cab9d09d8fc")},
                // See: GATKBAMIndex.getStartOfLastLinearBin(), BAMScheduler.advance(), IntervalOverlapFilteringIterator.advance()
                {new PRTest(b37KGReference, "unmappedFlagReadsInLastLinearBin.bam", "", "d7f23fd77d7dc7cb50d3397f644c6d8a")},
                {new PRTest(b37KGReference, "unmappedFlagReadsInLastLinearBin.bam", " -L 1", "c601db95b20248d012b0085347fcb6d1")},
                {new PRTest(b37KGReference, "unmappedFlagReadsInLastLinearBin.bam", " -L unmapped", "2d32440e47e8d9d329902fe573ad94ce")},
                {new PRTest(b37KGReference, "unmappedFlagReadsInLastLinearBin.bam", " -L 1 -L unmapped", "c601db95b20248d012b0085347fcb6d1")},
                {new PRTest(b37KGReference, "oneReadAllInsertion.bam", "",  "349650b6aa9e574b48a2a62627f37c7d")},
                {new PRTest(b37KGReference, "NA12878.1_10mb_2_10mb.bam", "",  "0c1cbe67296637a85e80e7a182f828ab")}
        };
    }

    @Test(dataProvider = "PRTest")
    public void testPrintReads(PRTest params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintReads" +
                        " -R " + params.reference +
                        " -I " + privateTestDir + params.bam +
                        params.args +
                        " --no_pg_tag" +
                        " -o %s",
                Arrays.asList(params.md5));
        executeTest("testPrintReads-"+params.args, spec).getFirst();
    }
}
