package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class MergeSegregatingAlternateAllelesIntegrationTest extends WalkerTest {

    public static String baseTestString(String reference, String VCF, int maxDist) {
        return "-T MergeSegregatingAlternateAlleles" +
                " -R " + reference +
                " --variant:vcf " + validationDataLocation + VCF +
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
                Arrays.asList("af5e1370822551c0c6f50f23447dc627"));
        executeTest("Merge sites within genomic distance of 1 [TEST ONE]", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 10)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("dd8c44ae1ef059a7fe85399467e102eb"));
        executeTest("Merge sites within genomic distance of 10 [TEST TWO]", spec);
    }

    @Test
    public void test3() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(hg18Reference, "merging_test_chr20_556259_756570.vcf", 100)
                        + " -L chr20:556259-756570",
                1,
                Arrays.asList("f81fd72ecaa57b3215406fcea860bcc5"));
        executeTest("Merge sites within genomic distance of 100 [TEST THREE]", spec);
    }


}