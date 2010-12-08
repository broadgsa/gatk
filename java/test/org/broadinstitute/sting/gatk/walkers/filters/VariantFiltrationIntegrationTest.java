package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class VariantFiltrationIntegrationTest extends WalkerTest {

    public static String baseTestString() {
        return "-T VariantFiltration -o %s -NO_HEADER -R " + b36KGReference;
    }


    @Test
    public void testNoAction() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B:variant,VCF " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("c9f49679f600c52e624c74cf65fd0a39"));
        executeTest("test no action", spec);
    }

    @Test
    public void testClusteredSnps() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -window 10 -B:variant,VCF " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("559021b518503df2ca1cc5cd18642248"));
        executeTest("test clustered SNPs", spec);
    }

    @Test
    public void testMask() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -mask foo -B:mask,VCF " + validationDataLocation + "vcfexample2.vcf -B:variant,VCF " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("467480f5c2e840b68d8ef591f1c562c7"));
        executeTest("test mask", spec);
    }

    @Test
    public void testFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo -B:variant,VCF " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("b48b5965d05f2d0f02c6f7261e67e9d2"));
        executeTest("test filter #1", spec);
    }

    @Test
    public void testFilter2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar -B:variant,VCF " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("42bfb7df7a83e5007a1d1da523a64d75"));
        executeTest("test filter #2", spec);
    }

    @Test
    public void testFilterWithSeparateNames() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterName ABF -filter 'AlleleBalance < 0.7' --filterName FSF -filter 'FisherStrand == 1.4' -B:variant,VCF " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("c2b1d53a50d99229b351470a39692eb6"));
        executeTest("test filter with separate names #2", spec);
    }

    @Test
    public void testGenotypeFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G_filter 'GQ == 0.60' -G_filterName foo -B:variant,VCF " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("6696e3f65a62ce912230d47cdb0c129b"));
        executeTest("test genotype filter #1", spec);
    }

    @Test
    public void testGenotypeFilter2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G_filter 'AF == 0.04 && isHomVar == 1' -G_filterName foo -B:variant,VCF " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("26e5b4ee954c9e0b5eb044afd4b88ee9"));
        executeTest("test genotype filter #2", spec);
    }
}
