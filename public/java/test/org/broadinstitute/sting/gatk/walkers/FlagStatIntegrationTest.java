package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class FlagStatIntegrationTest extends WalkerTest {

    @Test
    public void testFlagStat() {
        String md5 = "9c4039662f24bfd23ccf67973cb5df29";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T FlagStat -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000 -o %s",
                 1,
                 Arrays.asList(md5));
        executeTest("test flag stat", spec);
    }
}
