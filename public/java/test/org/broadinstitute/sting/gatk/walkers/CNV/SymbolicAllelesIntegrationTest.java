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
                " --no_cmdline_in_header";
    }


    @Test
    public void test1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(b36KGReference, "symbolic_alleles_1.vcf"),
                1,
                Arrays.asList("444a20659f67592a8284e0b7849e4302"));
        executeTest("Test symbolic alleles", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(b36KGReference, "symbolic_alleles_2.vcf"),
                1,
                Arrays.asList("93a24c019663a6011b4d6de12538df11"));
        executeTest("Test symbolic alleles mixed in with non-symbolic alleles", spec);
    }
}
