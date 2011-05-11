package org.broadinstitute.sting.playground.gatk.walkers.validation;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * tests for the ROD system in general; from rod system validation to empty VCF files
 */
public class RodSystemValidationIntegrationTest extends WalkerTest {

    public static String baseTestString1KG() {
            return "-T RodSystemValidation -o %s -R " + b36KGReference;
        }


    @Test
    public void testSimpleGeliPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " -B:eval,GeliText " + validationDataLocation + "ROD_validation/chr1.geli", 1,
                Arrays.asList("832efb29a6d4e8dbae374d3eeee17d9d"));
        executeTest("testSimpleGeliPileup", spec);
    }

    @Test
    public void testSimpleVCFPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " -B:eval,VCF3 " + validationDataLocation + "MultiSample.vcf", 1,
                Arrays.asList("7dc521914d1fd39cc0cab0220b6187f4"));
        executeTest("testSimpleVCFPileup", spec);
    }

    @Test
    public void testEmptyVCF() {
        File vcf = new File(validationDataLocation + "justHeader.vcf.idx");
        if (vcf.exists()) vcf.delete();

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " -B:eval,VCF3 " + validationDataLocation + "justHeader.vcf", 1,
                Arrays.asList("fa30619829498040414b45ec78151079"));
        executeTest("testEmptyVCF", spec);
    }


    @Test
    public void testComplexVCFPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " -B:eval,VCF3 " + validationDataLocation + "MultiSample.vcf" +
                " -B:eval2,VCF " + validationDataLocation + "NA12878.chr1_10mb_11mb.slx.indels.vcf4"
                , 1,
                Arrays.asList("ad9d697a015860089774452886cc1028"));
        executeTest("testComplexVCFPileup", spec);
    }

    @Test
    public void testBTIWithROD() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " -B:eval,VCF3 " + validationDataLocation + "MultiSample.vcf" +
                " -B:eval2,VCF " + validationDataLocation + "NA12878.chr1_10mb_11mb.slx.indels.vcf4" + " -BTI eval"
                , 1,
                Arrays.asList("92adff4129bc4c6abb96fb67b8bf61b6"));
        executeTest("testBTIWithROD", spec);
    }

    @Test
    public void testLargeComplexVCFPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " -B:eval,VCF3 " + validationDataLocation + "MultiSample.vcf" +
                " -B:eval2,VCF3 " + validationDataLocation + "CEU_hapmap_nogt_23.vcf" +
                " -B:eval3,VCF3 " + validationDataLocation + "CEU_hapmap_nogt_23.vcf" +
                " -L 1 -L 2 -L 20"
                , 1,
                Arrays.asList("2069cc9155af2ec10028acc8df4523b3"));
        executeTest("testLargeComplexVCFPileup", spec);
    }

    //@Test
    public void testBlockZippedVrsUnzippedVCF1() {
        final String vcfName = validationDataLocation + "bgzipped_vcfs/vcfexample.vcf";
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString1KG() + " -B:eval,VCF " + vcfName +
                " -B:eval2,VCF3 " + vcfName + ".gz" +
                " --PerLocusEqual"
                , 1,
                Arrays.asList("ab3da32eae65e8c15a9f4a787a190a37"));
        executeTest("testLargeComplexVCFPileup", spec);
    }
}
