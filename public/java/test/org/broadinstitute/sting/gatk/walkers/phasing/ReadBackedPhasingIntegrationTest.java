package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class ReadBackedPhasingIntegrationTest extends WalkerTest {

    public static String baseTestString(String reference, String reads, String VCF, int cacheWindowSize, int maxPhaseSites, double phaseQualityThresh) {
        return "-T ReadBackedPhasing" +
                " -R " + reference +
                " -I " + validationDataLocation + reads +
                " --variant " + validationDataLocation + VCF +
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
                Arrays.asList("0a41b96b04a87fdb99bc3342d48d2eba"));
        executeTest("MAX 10 het sites [TEST ONE]; require PQ >= 10", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 10, 10)
                        + " -L chr20:1232503-1332503",
                1,
                Arrays.asList("f7517896c899a872c24d8e823ac9deae"));
        executeTest("MAX 10 het sites [TEST TWO]; require PQ >= 10", spec);
    }

    @Test
    public void test3() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 2, 30)
                        + " -L chr20:332341-382503",
                1,
                Arrays.asList("cdbdd2f68c232012b6fe9a322b0ea24c"));
        executeTest("MAX 2 het sites [TEST THREE]; require PQ >= 30", spec);
    }

    @Test
    public void test4() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 5, 100)
                        + " -L chr20:332341-382503",
                1,
                Arrays.asList("6b70e3e4e28f9583d35d98bf8a7d0d59"));
        executeTest("MAX 5 het sites [TEST FOUR]; require PQ >= 100", spec);
    }

    @Test
    public void test5() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 1000, 7, 10)
                        + " -L chr20:332341-482503",
                1,
                Arrays.asList("6163a1fba27532da77765a7a11c55332"));
        executeTest("MAX 7 het sites [TEST FIVE]; require PQ >= 10; cacheWindow = 1000", spec);
    }

    @Test
    public void test6() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 10, 10)
                        + " -L chr20:652810-681757",
                1,
                Arrays.asList("61a7d05f9eb4317cf0e6937d72e1e7ec"));
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
