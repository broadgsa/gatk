package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
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
