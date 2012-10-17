package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

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



    private final static String SingleReadAligningOffChromosome1MD5 = "4a45fe1f85aaa8c4158782f2b6dee2bd";
    @Test
    public void testSingleReadAligningOffChromosome1() {
        String gatk_args = "-T Pileup "
                + " -I " + privateTestDir + "readOffb37contig1.bam"
                + " -R " + b37KGReference
                + " -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 1, Arrays.asList(SingleReadAligningOffChromosome1MD5));
        executeTest("Testing single read spanning off chromosome 1", spec);
    }

    @Test
    public void testSingleReadAligningOffChromosome1NoIndex() {
        String gatk_args = "-T Pileup "
                + " -I " + privateTestDir + "readOffb37contig1.noIndex.bam"
                + " -R " + b37KGReference
                + " -U ALLOW_UNINDEXED_BAM"
                + " -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 1, Arrays.asList(SingleReadAligningOffChromosome1MD5));
        executeTest("Testing single read spanning off chromosome 1 unindexed", spec);
    }
}
