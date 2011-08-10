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
                baseTestString() + " --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("8a105fa5eebdfffe7326bc5b3d8ffd1c"));
        executeTest("test no action", spec);
    }

    @Test
    public void testClusteredSnps() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -window 10 --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("27b13f179bb4920615dff3a32730d845"));
        executeTest("test clustered SNPs", spec);
    }

    @Test
    public void testMasks() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -maskName foo --mask:VCF3 " + validationDataLocation + "vcfexample2.vcf --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("578f9e774784c25871678e6464fd212b"));
        executeTest("test mask all", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -maskName foo --mask:VCF " + validationDataLocation + "vcfMask.vcf --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("bfa86a674aefca1b13d341cb14ab3c4f"));
        executeTest("test mask some", spec2);

        WalkerTestSpec spec3 = new WalkerTestSpec(
                baseTestString() + " -maskName foo -maskExtend 10 --mask:VCF " + validationDataLocation + "vcfMask.vcf --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("5939f80d14b32d88587373532d7b90e5"));
        executeTest("test mask extend", spec3);
    }

    @Test
    public void testFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("45219dbcfb6f81bba2ea0c35f5bfd368"));
        executeTest("test filter #1", spec);
    }

    @Test
    public void testFilter2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("c95845e817da7352b9b72bc9794f18fb"));
        executeTest("test filter #2", spec);
    }

    @Test
    public void testFilterWithSeparateNames() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterName ABF -filter 'AlleleBalance < 0.7' --filterName FSF -filter 'FisherStrand == 1.4' --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("b8cdd7f44ff1a395e0a9b06a87e1e530"));
        executeTest("test filter with separate names #2", spec);
    }

    @Test
    public void testGenotypeFilters() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'GQ == 0.60' -G_filterName foo --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("96b61e4543a73fe725e433f007260039"));
        executeTest("test genotype filter #1", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'AF == 0.04 && isHomVar == 1' -G_filterName foo --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("6c8112ab17ce39c8022c891ae73bf38e"));
        executeTest("test genotype filter #2", spec2);
    }

    @Test
    public void testDeletions() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterExpression 'QUAL < 100' --filterName foo --variant:VCF " + validationDataLocation + "twoDeletions.vcf", 1,
                Arrays.asList("569546fd798afa0e65c5b61b440d07ac"));
        executeTest("test deletions", spec);
    }
}
