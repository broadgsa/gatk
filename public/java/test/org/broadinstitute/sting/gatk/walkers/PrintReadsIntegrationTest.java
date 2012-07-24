package org.broadinstitute.sting.gatk.walkers;

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
                {new PRTest(hg18Reference, "HiSeq.1mb.bam", "", "dc8e5451dd29757c336013146010f73a")},
                {new PRTest(hg18Reference, "HiSeq.1mb.bam", " -compress 0", "fde82269c78c9e91e57286433531b4af")},
                {new PRTest(hg18Reference, "HiSeq.1mb.bam", " -simplifyBAM", "0531717b32a7e21c0de70b1526b0751f")},
                {new PRTest(hg18Reference, "HiSeq.1mb.bam", " -n 10", "cdc4ddf9ee1d2ecf37168da8ef23c270")},
                // See: GATKBAMIndex.getStartOfLastLinearBin(), BAMScheduler.advance(), IntervalOverlapFilteringIterator.advance()
                {new PRTest(b37KGReference, "unmappedFlagReadsInLastLinearBin.bam", "", "0a9ce949d07a84cb33a1a8e3358bf679")},
                {new PRTest(b37KGReference, "unmappedFlagReadsInLastLinearBin.bam", " -L 1", "6e920b8505e7e95d67634b0905237dbc")},
                {new PRTest(b37KGReference, "unmappedFlagReadsInLastLinearBin.bam", " -L unmapped", "13bb9a91b1d4dd2425f73302b8a1ac1c")},
                {new PRTest(b37KGReference, "unmappedFlagReadsInLastLinearBin.bam", " -L 1 -L unmapped", "6e920b8505e7e95d67634b0905237dbc")},
                {new PRTest(b37KGReference, "oneReadAllInsertion.bam", "",  "6caec4f8a25befb6aba562955401af93")}
        };
    }

    @Test(dataProvider = "PRTest")
    public void testPrintReads(PRTest params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintReads" +
                        " -R " + params.reference +
                        " -I " + privateTestDir + params.bam +
                        params.args +
                        " -o %s",
                Arrays.asList(params.md5));
        executeTest("testPrintReads-"+params.args, spec).getFirst();
    }
}
