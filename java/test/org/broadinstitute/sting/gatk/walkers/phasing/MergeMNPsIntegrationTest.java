package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class MergeMNPsIntegrationTest extends WalkerTest {

    public static String baseTestString(String reference, String VCF, int maxDistMNP) {
        return "-T MergeMNPs" +
                " -R " + reference +
                " -B:variant,VCF " + validationDataLocation + VCF +
                " --maxGenomicDistanceForMNP " + maxDistMNP +
                " -o %s" +
                " -NO_HEADER";
    }


    @Test
    public void test1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 1)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("e312b7d3854d5b2834a370659514a813"));
        executeTest("Merge MNP sites within genomic distance of 1 [TEST ONE]", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 10)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("681f50e45f1d697370d2c355df2e18bc"));
        executeTest("Merge MNP sites within genomic distance of 10 [TEST TWO]", spec);
    }

    @Test
    public void test3() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 100)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("0bccb0ef928a108418246bec01098083"));
        executeTest("Merge MNP sites within genomic distance of 100 [TEST THREE]", spec);
    }


}