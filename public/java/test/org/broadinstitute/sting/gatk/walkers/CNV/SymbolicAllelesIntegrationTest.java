package org.broadinstitute.sting.gatk.walkers.CNV;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class SymbolicAllelesIntegrationTest extends WalkerTest {

    public static String baseTestString(String reference, String VCF) {
        return "-T CombineVariants" +
                " -R " + reference +
                " --variant:vcf " + testDir + VCF +
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
                Arrays.asList("c79137da24ad4dc15cedc742de39247f"));
        executeTest("Test symbolic alleles", spec);
    }

    @Test
    public void test2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(b36KGReference, "symbolic_alleles_2.vcf"),
                1,
                Arrays.asList("3f6cbbd5fdf164d87081a3af19eeeba7"));
        executeTest("Test symbolic alleles mixed in with non-symbolic alleles", spec);
    }
}
