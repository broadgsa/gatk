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
                baseTestString() + " -B:variant,VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("4cc077eb3d343e6b7ba12bff86ebe347"));
        executeTest("test no action", spec);
    }

    @Test
    public void testClusteredSnps() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -window 10 -B:variant,VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("ada5540bb3d9b6eb8f1337ba01e90a94"));
        executeTest("test clustered SNPs", spec);
    }

    @Test
    public void testMasks() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -mask foo -B:mask,VCF3 " + validationDataLocation + "vcfexample2.vcf -B:variant,VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("b0fcac4af3526e3b2a37602ab4c0e6ae"));
        executeTest("test mask all", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -mask foo -B:mask,VCF " + validationDataLocation + "vcfMask.vcf -B:variant,VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("b64baabe905a5d197cc1ab594147d3d5"));
        executeTest("test mask some", spec2);

        WalkerTestSpec spec3 = new WalkerTestSpec(
                baseTestString() + " -mask foo -maskExtend 10 -B:mask,VCF " + validationDataLocation + "vcfMask.vcf -B:variant,VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("0eff92fe72024d535c44b98e1e9e1993"));
        executeTest("test mask extend", spec3);
    }

    @Test
    public void testFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo -B:variant,VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("7a40795147cbfa92941489d7239aad92"));
        executeTest("test filter #1", spec);
    }

    @Test
    public void testFilter2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar -B:variant,VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("e9dd4991b1e325847c77d053dfe8ee54"));
        executeTest("test filter #2", spec);
    }

    @Test
    public void testFilterWithSeparateNames() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterName ABF -filter 'AlleleBalance < 0.7' --filterName FSF -filter 'FisherStrand == 1.4' -B:variant,VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("9ded2cce63b8d97550079047051d80a3"));
        executeTest("test filter with separate names #2", spec);
    }

    @Test
    public void testGenotypeFilters() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'GQ == 0.60' -G_filterName foo -B:variant,VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("6696e3f65a62ce912230d47cdb0c129b"));
        executeTest("test genotype filter #1", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'AF == 0.04 && isHomVar == 1' -G_filterName foo -B:variant,VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("26e5b4ee954c9e0b5eb044afd4b88ee9"));
        executeTest("test genotype filter #2", spec2);
    }

    @Test
    public void testDeletions() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterExpression 'QUAL < 100' --filterName foo -B:variant,VCF " + validationDataLocation + "twoDeletions.vcf", 1,
                Arrays.asList("e63b58be33c9126ad6cc55489aac539b"));
        executeTest("test deletions", spec);
    }
}
