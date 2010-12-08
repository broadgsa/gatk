package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class MergeSegregatingAlternateAllelesIntegrationTest extends WalkerTest {

    public static String baseTestString(String reference, String VCF, int maxDist) {
        return "-T MergeSegregatingAlternateAlleles" +
                " -R " + reference +
                " -B:variant,VCF " + validationDataLocation + VCF +
                " --maxGenomicDistance " + maxDist +
                " -o %s" +
                " -NO_HEADER";
    }


    @Test
    public void test1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 1)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("e16f957d888054ae0518e25660295241"));
        executeTest("Merge sites within genomic distance of 1 [TEST ONE]", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 10)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("122a482090677c7619c2105d44e00d11"));
        executeTest("Merge sites within genomic distance of 10 [TEST TWO]", spec);
    }

    @Test
    public void test3() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 100)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("bc6a8c8a42bb2601db98e88e9ad74748"));
        executeTest("Merge sites within genomic distance of 100 [TEST THREE]", spec);
    }


}