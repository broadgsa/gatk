package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class MergeSegregatingPolymorphismsIntegrationTest extends WalkerTest {

    public static String baseTestString(String reference, String VCF, int maxDistMNP) {
        return "-T MergeSegregatingPolymorphisms" +
                " -R " + reference +
                " -B:variant,VCF " + validationDataLocation + VCF +
                " --maxGenomicDistanceForMNP " + maxDistMNP +
                " -o %s";
    }


    @Test
    public void test1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 1)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("38d88fc6c1880e76ce402cfc60669726"));
        executeTest("Merge MNP het sites within genomic distance of 1 [TEST ONE]", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 10)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("e14a2f062391b6c3f8b36b1b4eed628b"));
        executeTest("Merge MNP het sites within genomic distance of 10 [TEST TWO]", spec);
    }

    @Test
    public void test3() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 100)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("fdb4ad8ced0cb461deaff1555258008e"));
        executeTest("Merge MNP het sites within genomic distance of 100 [TEST THREE]", spec);
    }


}