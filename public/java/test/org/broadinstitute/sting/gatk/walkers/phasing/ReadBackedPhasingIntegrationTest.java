package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class ReadBackedPhasingIntegrationTest extends WalkerTest {

    public static String baseTestString(String reference, String reads, String VCF, int cacheWindowSize, int maxPhaseSites, double phaseQualityThresh) {
        return "-T ReadBackedPhasing" +
                " -R " + reference +
                " -I " + validationDataLocation + reads +
                " -B:variant,VCF " + validationDataLocation + VCF +
                " --cacheWindowSize " + cacheWindowSize +
                " --maxPhaseSites " + maxPhaseSites +
                " --phaseQualityThresh " + phaseQualityThresh +
                " -o %s" +
                " -NO_HEADER";
    }


    @Test
    public void test1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 10, 10)
                        + " -L chr20:332341-382503",
                1,
                Arrays.asList("6020a68bbec97fcd87819c10cd4e2470"));
        executeTest("MAX 10 het sites [TEST ONE]; require PQ >= 10", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 10, 10)
                        + " -L chr20:1232503-1332503",
                1,
                Arrays.asList("712c2145df4756c9a15758865d8007b5"));
        executeTest("MAX 10 het sites [TEST TWO]; require PQ >= 10", spec);
    }

    @Test
    public void test3() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 2, 30)
                        + " -L chr20:332341-382503",
                1,
                Arrays.asList("297e0896e4761529d979f40f5ad694db"));
        executeTest("MAX 2 het sites [TEST THREE]; require PQ >= 30", spec);
    }

    @Test
    public void test4() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 5, 100)
                        + " -L chr20:332341-382503",
                1,
                Arrays.asList("52a17f14692d726d3b726cf0ae7f2a09"));
        executeTest("MAX 5 het sites [TEST FOUR]; require PQ >= 100", spec);
    }

    @Test
    public void test5() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 1000, 7, 10)
                        + " -L chr20:332341-482503",
                1,
                Arrays.asList("af768f7958b8f4599c2374f1cc2fc613"));
        executeTest("MAX 7 het sites [TEST FIVE]; require PQ >= 10; cacheWindow = 1000", spec);
    }

    @Test
    public void test6() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "phasing_test_chr20_332341_1332503.bam", "phasing_test_chr20_332341_1332503.vcf", 20000, 10, 10)
                        + " -L chr20:652810-681757",
                1,
                Arrays.asList("3dd886672f59a47908b94136d0427bb0"));
        executeTest("MAX 10 het sites [TEST SIX]; require PQ >= 10; cacheWindow = 20000; has inconsistent sites", spec);
    }

}
