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
                Arrays.asList("ee5e4e00bf25a912e8cab3e768fa0e7d"));
        executeTest("test clustered SNPs", spec);
    }

    @Test
    public void testMask() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -mask foo -B mask,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("77451c3ba8a070343e69157cdaf2be92"));
        executeTest("test mask", spec);
    }

    @Test
    public void testFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("3b44be79d676536bd6c0f32774091fee"));
        executeTest("test filter #1", spec);
    }

    @Test
    public void testFilter2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("baa425fa8e0761a05733a1a9d62d02ff"));
        executeTest("test filter #2", spec);
    }
}