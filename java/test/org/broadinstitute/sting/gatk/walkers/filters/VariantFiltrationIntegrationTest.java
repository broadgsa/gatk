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
                Arrays.asList("e8e1898d65eb77ec32b7eca6764864f3"));
        executeTest("test no action", spec);
    }

    @Test
    public void testClusteredSnps() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -window 10 -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("25da669bc6411d31f585b01be85e8429"));
        executeTest("test clustered SNPs", spec);
    }

    @Test
    public void testMask() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -mask foo -B mask,VCF," + validationDataLocation + "vcfexample2.vcf -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("bd08169351c4cc9a2d95f33a2a73c7fa"));
        executeTest("test mask", spec);
    }

    @Test
    public void testFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("9e7631aeda1c174f5aac6712279ed861"));
        executeTest("test filter #1", spec);
    }

    @Test
    public void testFilter2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("9009742f5697c5a85f9e41e1d7d6f260"));
        executeTest("test filter #2", spec);
    }

    @Test
    public void testFilterWithSeparateNames() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterName ABF -filter 'AlleleBalance < 70.0' --filterName FSF -filter 'FisherStrand == 1.4' -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("6b280c5abd6e911df8e4fca58c7b216b"));
        executeTest("test filter with separate names #2", spec);
    }
}