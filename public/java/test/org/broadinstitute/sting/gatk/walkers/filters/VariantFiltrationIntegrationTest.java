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
                baseTestString() + " --variant " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("49471b44ac165929d3ff81f98ce19063"));
        executeTest("test no action", spec);
    }

    @Test
    public void testClusteredSnps() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -window 10 --variant " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("8b45895d7ae1f36b70e7fd26aa9451d3"));
        executeTest("test clustered SNPs", spec);
    }

    @Test
    public void testMask1() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -maskName foo --mask " + testDir + "vcfexample2.vcf --variant " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("06307029f5da87ae4edd9804063a98f9"));
        executeTest("test mask all", spec1);
    }

    @Test
    public void testMask2() {
        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -maskName foo --mask:VCF " + testDir + "vcfMask.vcf --variant " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("1fd06f6b2642685093ed36342f002b58"));
        executeTest("test mask some", spec2);
    }

    @Test
    public void testMask3() {
        WalkerTestSpec spec3 = new WalkerTestSpec(
                baseTestString() + " -maskName foo -maskExtend 10 --mask:VCF " + testDir + "vcfMask.vcf --variant " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("d8c5206d5d13477a5929fb1ae5a6bfc4"));
        executeTest("test mask extend", spec3);
    }

    @Test
    public void testFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo --variant " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("a3be095e8aa75d9ef4235b9487527307"));
        executeTest("test filter #1", spec);
    }

    @Test
    public void testFilter2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar --variant " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("bd1361ddc52d73b8cd7adeb9e5c47200"));
        executeTest("test filter #2", spec);
    }

    @Test
    public void testFilterWithSeparateNames() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterName ABF -filter 'AlleleBalance < 0.7' --filterName FSF -filter 'FisherStrand == 1.4' --variant " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("4a43ec0285433df426ab482f88cf7ca6"));
        executeTest("test filter with separate names #2", spec);
    }

    @Test
    public void testGenotypeFilters1() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'GQ == 0.60' -G_filterName foo --variant " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("5ee4485a022e163645c08b9691384f67"));
        executeTest("test genotype filter #1", spec1);
    }

    @Test
    public void testGenotypeFilters2() {
        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'AF == 0.04 && isHomVar == 1' -G_filterName foo --variant " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("d0a068c8cfb0758d2a8d471383f39b68"));
        executeTest("test genotype filter #2", spec2);
    }

    @Test
    public void testDeletions() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterExpression 'QUAL < 100' --filterName foo --variant:VCF " + testDir + "twoDeletions.vcf", 1,
                Arrays.asList("a1c02a5a90f1262e9eb3d2cad1fd08f2"));
        executeTest("test deletions", spec);
    }
}
