package org.broadinstitute.sting.playground.gatk.walkers.reducereads;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public class ReduceReadsIntegrationTest extends WalkerTest {
    final static String REF = b37KGReference;
    final String BAM = validationDataLocation + "NA12878.HiSeq.b37.chr20.10_11mb.bam";
    final String L = "20:10,100,000-10,200,000";

    private void RRTest(String args, String md5) {
        String base = String.format("-T ReduceReads -R %s -I %s -L %s", REF, BAM, L) + " -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(base + args, Arrays.asList(md5));
        executeTest("testReduceReads1: args=" + args, spec);
    }

    @Test()
    public void testReduceReadsBase() {
        RRTest("", "9aac475601d80eeb02400cbc73272b26");
    }

    @Test()
    public void testReduceReads50MaxReads() {
        RRTest(" -mravs 50", "eb2d8c2f1e66d7d0bf767ac55420027e");
    }

    @Test()
    public void testReduceReadsMinBasesForConsensus10000() {
        RRTest(" -mbrc 10000", "b42706d9a2621b9b63502704af00e0da");
    }

}

