package org.broadinstitute.sting.gatk.walkers.CNV;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class SymbolicAllelesIntegrationTest extends WalkerTest {

    public static String baseTestString(String reference, String VCF) {
        return "-T CombineVariants" +
                " -R " + reference +
                " --variant:vcf " + validationDataLocation + VCF +
                " -filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED" +
                " -genotypeMergeOptions REQUIRE_UNIQUE" +
                " -setKey null" +
                " -o %s" +
                " -NO_HEADER";
    }


    @Test
    public void test1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(b36KGReference, "symbolic_alleles_1.vcf"),
                1,
                Arrays.asList("8c44f58039ec444bc376fa732d4856dd"));
        executeTest("Test symbolic alleles", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(b36KGReference, "symbolic_alleles_2.vcf"),
                1,
                Arrays.asList("d2fc144d5996530b9161af3d22af41d9"));
        executeTest("Test symbolic alleles mixed in with non-symbolic alleles", spec);
    }
}
