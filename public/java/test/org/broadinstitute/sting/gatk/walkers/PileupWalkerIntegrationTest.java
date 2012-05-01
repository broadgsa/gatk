package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Dec 1, 2009
 * Time: 9:03:34 AM
 * To change this template use File | Settings | File Templates.
 */
public class PileupWalkerIntegrationTest extends WalkerTest {

    @Test
    public void testGnarleyFHSPileup() {
        String gatk_args = "-T Pileup -I " + validationDataLocation + "FHS_Pileup_Test.bam "
                 + "-R " + hg18Reference
                 +  " -L chr15:46,347,148 -o %s";
        String expected_md5 = "526c93b0fa660d6b953b57103e59fe98";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 1, Arrays.asList(expected_md5));
        executeTest("Testing the standard (no-indel) pileup on three merged FHS pools with 27 deletions in 969 bases", spec);
    }
}
