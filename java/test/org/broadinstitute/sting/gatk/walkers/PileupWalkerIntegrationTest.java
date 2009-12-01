package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

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
        String gatk_args = "-T Pileup -I /humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_Pileup_Test.bam "
                 + "-R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"
                 +  " -L chr15:46,347,148 -o %s";
        String expected_md5 = "98040c47b623bb133cd296ce12768d49";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 1, Arrays.asList(expected_md5));
        executeTest("Testing the standard (no-indel) pileup on three merged FHS pools with 27 deletions in 969 bases", spec);
    }
}
