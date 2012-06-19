package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class ReadBackedPhasingIntegrationTest extends WalkerTest {

    public static String baseTestString(String reference, String reads, String VCF, int cacheWindowSize, int maxPhaseSites, double phaseQualityThresh) {
        return "-T ReadBackedPhasing" +
                " -R " + reference +
                " -I " + validationDataLocation + reads +
                " --variant " + ( VCF.contains("phasing_test") ? testDir : validationDataLocation) + VCF +
                " --cacheWindowSize " + cacheWindowSize +
                " --maxPhaseSites " + maxPhaseSites +
                " --phaseQualityThresh " + phaseQualityThresh +
                " -o %s" +
                " --no_cmdline_in_header";
    }


    @Test
    public void test1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 10, 10)
                        + " -L chr20:332341-382503",
                1,
                Arrays.asList("442c819569417c1b7d6be9f41ce05394"));
        executeTest("MAX 10 het sites [TEST ONE]; require PQ >= 10", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 10, 10)
                        + " -L chr20:1232503-1332503",
                1,
                Arrays.asList("2a51ee7d3c024f2410dcee40c5412993"));
        executeTest("MAX 10 het sites [TEST TWO]; require PQ >= 10", spec);
    }

    @Test
    public void test3() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 2, 30)
                        + " -L chr20:332341-382503",
                1,
                Arrays.asList("85bc9b03e24159f746dbd0cb988f9ec8"));
        executeTest("MAX 2 het sites [TEST THREE]; require PQ >= 30", spec);
    }

    @Test
    public void test4() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 5, 100)
                        + " -L chr20:332341-382503",
                1,
                Arrays.asList("96bb413a83c777ebbe622438e4565e8f"));
        executeTest("MAX 5 het sites [TEST FOUR]; require PQ >= 100", spec);
    }

    @Test
    public void test5() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 1000, 7, 10)
                        + " -L chr20:332341-482503",
                1,
                Arrays.asList("7d2402f055d243e2208db9ea47973e13"));
        executeTest("MAX 7 het sites [TEST FIVE]; require PQ >= 10; cacheWindow = 1000", spec);
    }

    @Test
    public void test6() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 10, 10)
                        + " -L chr20:652810-681757",
                1,
                Arrays.asList("72682b3f27c33580d2d4515653ba6de7"));
        executeTest("MAX 10 het sites [TEST SIX]; require PQ >= 10; cacheWindow = 20000; has inconsistent sites", spec);
    }

    @Test
    public void test7() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "CEU.trio.2010_03.genotypes.hg18.vcf", 20000, 10, 10)
                        + " -L chr20:332341-802503",
                1,
                Arrays.asList("44eb225ab3167651ec0a9e1fdcc83d34"));
        executeTest("Use trio-phased VCF, but ignore its phasing [TEST SEVEN]", spec);
    }

    @Test
    public void test8() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "CEU.trio.2010_03.genotypes.hg18.vcf", 20000, 10, 10)
                        + " -L chr20:332341-802503" + " -respectPhaseInInput",
                1,
                Arrays.asList("e3549b89d49092e73cc6eb21f233471c"));
        executeTest("Use trio-phased VCF, and respect its phasing [TEST EIGHT]", spec);
    }

}
