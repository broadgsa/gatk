package org.broadinstitute.sting.playground.gatk.walkers.validation;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

/**
 * The pile-up tests, that test any changes to the underlying ROD system
 */
public class RodSystemValidationIntegrationTest extends WalkerTest {

    public static String baseTestString1KG() {
            return "-T RodSystemValidation -o %s -R " + oneKGLocation + "reference/human_b36_both.fasta";
        }


    @Test
    public void testSimpleGeliPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " -B eval,GeliText," + validationDataLocation + "ROD_validation/chr1.geli", 1,
                Arrays.asList("832efb29a6d4e8dbae374d3eeee17d9d"));
        executeTest("testSimpleGeliPileup", spec);
    }

    @Test
    public void testSimpleVCFPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " -B eval,VCF," + validationDataLocation + "MultiSample.vcf", 1,
                Arrays.asList("5d84c75746738833b6c9441d9d614553"));
        executeTest("testSimpleVCFPileup", spec);
    }

    @Test
    public void testComplexVCFPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " -B eval,VCF," + validationDataLocation + "MultiSample.vcf" +
                " -B eval,VCF," + validationDataLocation + "NA12878.chr1_10mb_11mb.slx.indels.vcf"
                , 1,
                Arrays.asList("6dd0ed0a6fe7096ccb66beffb8d455da"));
        executeTest("testComplexVCFPileup", spec);
    }

    @Test
    public void testLargeComplexVCFPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " -B eval,VCF," + validationDataLocation + "MultiSample.vcf" +
                " -B eval,VCF," + validationDataLocation + "NA12878.chr1_10mb_11mb.slx.indels.vcf" +
                " -B eval,VCF," + validationDataLocation + "CEU_hapmap_nogt_23.vcf" +
                " -L 1 -L 2 -L 20"
                , 1,
                Arrays.asList("ab3da32eae65e8c15a9f4a787a190a37"));
        executeTest("testLargeComplexVCFPileup", spec);
    }
}
