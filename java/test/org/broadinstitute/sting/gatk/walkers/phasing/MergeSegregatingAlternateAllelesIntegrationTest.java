package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class MergeSegregatingAlternateAllelesIntegrationTest extends WalkerTest {

    public static String baseTestString(String reference, String VCF, int maxDistMNP) {
        return "-T MergeSegregatingAlternateAlleles" +
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
                Arrays.asList("e6a14fc97dbd0aaa8e6a4d9a7f1616a6"));
        executeTest("Merge MNP het sites within genomic distance of 1 [TEST ONE]", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 10)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("cc2b45c85a51b4998e30758c48f61940"));
        executeTest("Merge MNP het sites within genomic distance of 10 [TEST TWO]", spec);
    }

    @Test
    public void test3() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 100)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("47300cc7a5a7d84b3c279f04c4567739"));
        executeTest("Merge MNP het sites within genomic distance of 100 [TEST THREE]", spec);
    }


}