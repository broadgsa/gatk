package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class VariantFiltrationIntegrationTest extends WalkerTest {

    public static String baseTestString() {
        return "-T VariantFiltration -o %s -R /broad/1KG/reference/human_b36_both.fasta";
    }


    @Test
    public void testNoAction() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("49817e684effce7b6f3d5776dc781988"));
        executeTest("test no action", spec);
    }

    @Test
    public void testClusteredSnps() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -window 10 -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("7a1e834a372c3f3acfde82942be8c86c"));
        executeTest("test clustered SNPs", spec);
    }

    @Test
    public void testMask() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -mask foo -B mask,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("6b33aab1f199bd7856256dabdf729e1a"));
        executeTest("test mask", spec);
    }

    @Test
    public void testFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("7757178f6337e4c49097f62e8949f986"));
        executeTest("test filter #1", spec);
    }

    @Test
    public void testFilter2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("25354adffa90857be2900a0369cb2709"));
        executeTest("test filter #2", spec);
    }
}