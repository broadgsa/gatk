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
                Arrays.asList("f0e9796bc9866201aed0e97d76ed4a84"));
        executeTest("Merge MNP het sites within genomic distance of 1 [TEST ONE]", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 10)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("ab2b05dc6e7e4464e8f9e08d1d02f8ae"));
        executeTest("Merge MNP het sites within genomic distance of 10 [TEST TWO]", spec);
    }

    @Test
    public void test3() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 100)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("d446a0076f9d3d13c128bc8402b087a4"));
        executeTest("Merge MNP het sites within genomic distance of 100 [TEST THREE]", spec);
    }


}