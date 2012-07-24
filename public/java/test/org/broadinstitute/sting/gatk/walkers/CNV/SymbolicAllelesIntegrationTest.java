package org.broadinstitute.sting.gatk.walkers.CNV;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class SymbolicAllelesIntegrationTest extends WalkerTest {

    public static String baseTestString(String reference, String VCF) {
        return "-T CombineVariants" +
                " -R " + reference +
                " --variant:vcf " + privateTestDir + VCF +
                " -filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED" +
                " -genotypeMergeOptions REQUIRE_UNIQUE" +
                " -setKey null" +
                " -o %s" +
                " --no_cmdline_in_header";
    }

    @Test(enabled = true)
    public void test1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(b36KGReference, "symbolic_alleles_1.vcf"),
                1,
                Arrays.asList("5bafc5a99ea839e686e55de93f91fd5c"));
        executeTest("Test symbolic alleles", spec);
    }

    @Test(enabled = true)
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(b36KGReference, "symbolic_alleles_2.vcf"),
                1,
                Arrays.asList("bf5a09f783ab1fa44774c81f91d10921"));
        executeTest("Test symbolic alleles mixed in with non-symbolic alleles", spec);
    }
}
