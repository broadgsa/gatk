package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class VariantFiltrationIntegrationTest extends WalkerTest {

    public static String baseTestString() {
        return "-T VariantFiltration -o %s -R " + oneKGLocation + "reference/human_b36_both.fasta";
    }


    @Test
    public void testNoAction() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("e0543c72ed36f4c0c43d791ad44aa96a"));
        executeTest("test no action", spec);
    }

    @Test
    public void testClusteredSnps() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -window 10 -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("2385975931cd06fca452655bebf5c379"));
        executeTest("test clustered SNPs", spec);
    }

    @Test
    public void testMask() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -mask foo -B mask,VCF," + validationDataLocation + "vcfexample2.vcf -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("6743efa09985206819adcf8eaf5ff936"));
        executeTest("test mask", spec);
    }

    @Test
    public void testFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("e054f57d3794ce8a57cae92f16886cf0"));
        executeTest("test filter #1", spec);
    }

    @Test
    public void testFilter2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("db09a1f7ff523087ea0f6a56f3febfe7"));
        executeTest("test filter #2", spec);
    }

    @Test
    public void testFilterWithSeparateNames() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterName ABF -filter 'AlleleBalance < 70.0' --filterName FSF -filter 'FisherStrand == 1.4' -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("edbd505d8d55b4ba71f99d7006871db0"));
        executeTest("test filter with separate names #2", spec);
    }
}