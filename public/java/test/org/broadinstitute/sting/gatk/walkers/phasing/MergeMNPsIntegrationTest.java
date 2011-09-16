package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class MergeMNPsIntegrationTest extends WalkerTest {

    public static String baseTestString(String reference, String VCF, int maxDistMNP) {
        return "-T MergeMNPs" +
                " -R " + reference +
                " --variant:vcf " + validationDataLocation + VCF +
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
                Arrays.asList("7f11f7f75d1526077f0173c7ed1fc6c4"));
        executeTest("Merge MNP sites within genomic distance of 1 [TEST ONE]", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 10)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("53dd312468296826bdd3c22387390c88"));
        executeTest("Merge MNP sites within genomic distance of 10 [TEST TWO]", spec);
    }

    @Test
    public void test3() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 100)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("e26f92d2fb9f4eaeac7f9d8ee27410ee"));
        executeTest("Merge MNP sites within genomic distance of 100 [TEST THREE]", spec);
    }


}