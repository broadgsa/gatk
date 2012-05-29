package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class VariantFiltrationIntegrationTest extends WalkerTest {

    public static String baseTestString() {
        return "-T VariantFiltration -o %s --no_cmdline_in_header -R " + b36KGReference;
    }


    @Test
    public void testNoAction() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("5720826c2bf6cbc762e4a888ef58c3f2"));
        executeTest("test no action", spec);
    }

    @Test
    public void testClusteredSnps() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -window 10 --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("d7c2a4b0c1b2b982847508997ba57ebf"));
        executeTest("test clustered SNPs", spec);
    }

    @Test
    public void testMask1() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -maskName foo --mask:VCF3 " + testDir + "vcfexample2.vcf --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("890774962576c407d8a17ed57cf704c1"));
        executeTest("test mask all", spec1);
    }

    @Test
    public void testMask2() {
        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -maskName foo --mask:VCF " + testDir + "vcfMask.vcf --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("8864573dbf52908501140e6b0afcbc90"));
        executeTest("test mask some", spec2);
    }

    @Test
    public void testMask3() {
        WalkerTestSpec spec3 = new WalkerTestSpec(
                baseTestString() + " -maskName foo -maskExtend 10 --mask:VCF " + testDir + "vcfMask.vcf --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("42a1c08763f151073a49e3c7bb68028b"));
        executeTest("test mask extend", spec3);
    }

    @Test
    public void testFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("ef8100c3b7c67d28571cbda771c414c2"));
        executeTest("test filter #1", spec);
    }

    @Test
    public void testFilter2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("318ed3874fd42b7da8c59554a25a1fab"));
        executeTest("test filter #2", spec);
    }

    @Test
    public void testFilterWithSeparateNames() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterName ABF -filter 'AlleleBalance < 0.7' --filterName FSF -filter 'FisherStrand == 1.4' --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("9cb398e78a38a7bc5e839e28c8dae2eb"));
        executeTest("test filter with separate names #2", spec);
    }

    @Test
    public void testGenotypeFilters1() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'GQ == 0.60' -G_filterName foo --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("b38709f932b969e4267603333863269e"));
        executeTest("test genotype filter #1", spec1);
    }

    @Test
    public void testGenotypeFilters2() {
        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'AF == 0.04 && isHomVar == 1' -G_filterName foo --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("0e1457e678326e44e92ee13e84414e0f"));
        executeTest("test genotype filter #2", spec2);
    }

    @Test
    public void testDeletions() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterExpression 'QUAL < 100' --filterName foo --variant:VCF " + testDir + "twoDeletions.vcf", 1,
                Arrays.asList("569546fd798afa0e65c5b61b440d07ac"));
        executeTest("test deletions", spec);
    }
}
