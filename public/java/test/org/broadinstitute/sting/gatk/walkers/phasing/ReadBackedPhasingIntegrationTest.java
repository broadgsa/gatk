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
                Arrays.asList("2520f93505fda28d44f618a0123d593b"));
        executeTest("MAX 10 het sites [TEST ONE]; require PQ >= 10", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 10, 10)
                        + " -L chr20:1232503-1332503",
                1,
                Arrays.asList("965b8f448365b7f4a124d32e809eb048"));
        executeTest("MAX 10 het sites [TEST TWO]; require PQ >= 10", spec);
    }

    @Test
    public void test3() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 2, 30)
                        + " -L chr20:332341-382503",
                1,
                Arrays.asList("60f5bb699335f47cdc505322c5be3803"));
        executeTest("MAX 2 het sites [TEST THREE]; require PQ >= 30", spec);
    }

    @Test
    public void test4() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 5, 100)
                        + " -L chr20:332341-382503",
                1,
                Arrays.asList("023c2fb43b50807cfd46841ed6f0d215"));
        executeTest("MAX 5 het sites [TEST FOUR]; require PQ >= 100", spec);
    }

    @Test
    public void test5() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 1000, 7, 10)
                        + " -L chr20:332341-482503",
                1,
                Arrays.asList("e5e6e9f84d108d5b001aa53017d2801e"));
        executeTest("MAX 7 het sites [TEST FIVE]; require PQ >= 10; cacheWindow = 1000", spec);
    }

    @Test
    public void test6() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 10, 10)
                        + " -L chr20:652810-681757",
                1,
                Arrays.asList("8fc53bfbea2754ff8577460786a3400c"));
        executeTest("MAX 10 het sites [TEST SIX]; require PQ >= 10; cacheWindow = 20000; has inconsistent sites", spec);
    }

    @Test
    public void test7() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "CEU.trio.2010_03.genotypes.hg18.vcf", 20000, 10, 10)
                        + " -L chr20:332341-802503",
                1,
                Arrays.asList("c37548b333b65f58d0edfc5c2a62a28a"));
        executeTest("Use trio-phased VCF, but ignore its phasing [TEST SEVEN]", spec);
    }

    @Test
    public void test8() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "CEU.trio.2010_03.genotypes.hg18.vcf", 20000, 10, 10)
                        + " -L chr20:332341-802503" + " -respectPhaseInInput",
                1,
                Arrays.asList("dfc7cdddd702e63d46d04f61a3ecd720"));
        executeTest("Use trio-phased VCF, and respect its phasing [TEST EIGHT]", spec);
    }

}
