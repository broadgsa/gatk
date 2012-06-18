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
                Arrays.asList("dfa5dff09fa964b06da19c0f4aff6928"));
        executeTest("test no action", spec);
    }

    @Test
    public void testClusteredSnps() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -window 10 --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("4a4596929f9fe983d8868ca142567781"));
        executeTest("test clustered SNPs", spec);
    }

    @Test
    public void testMask1() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -maskName foo --mask:VCF3 " + testDir + "vcfexample2.vcf --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("1719462cd17986c33e59e45b69df0270"));
        executeTest("test mask all", spec1);
    }

    @Test
    public void testMask2() {
        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -maskName foo --mask:VCF " + testDir + "vcfMask.vcf --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("db19ff7d90c82cda09fb3c3878100eb5"));
        executeTest("test mask some", spec2);
    }

    @Test
    public void testMask3() {
        WalkerTestSpec spec3 = new WalkerTestSpec(
                baseTestString() + " -maskName foo -maskExtend 10 --mask:VCF " + testDir + "vcfMask.vcf --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("a9e417cba21585c786d4b9930265ea31"));
        executeTest("test mask extend", spec3);
    }

    @Test
    public void testFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("4160904b180d1f62a6bf50de6728ce00"));
        executeTest("test filter #1", spec);
    }

    @Test
    public void testFilter2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("df80db30c7836731ac7c8c3d4fc005b4"));
        executeTest("test filter #2", spec);
    }

    @Test
    public void testFilterWithSeparateNames() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterName ABF -filter 'AlleleBalance < 0.7' --filterName FSF -filter 'FisherStrand == 1.4' --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("71ce6c0952831cb68f575aa0173dce2b"));
        executeTest("test filter with separate names #2", spec);
    }

    @Test
    public void testGenotypeFilters1() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'GQ == 0.60' -G_filterName foo --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("179f7f2a90c0e6c656109aac9b775476"));
        executeTest("test genotype filter #1", spec1);
    }

    @Test
    public void testGenotypeFilters2() {
        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'AF == 0.04 && isHomVar == 1' -G_filterName foo --variant:VCF3 " + testDir + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("22e07c27feb9017a130dfb045c5b29b9"));
        executeTest("test genotype filter #2", spec2);
    }

    @Test
    public void testDeletions() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterExpression 'QUAL < 100' --filterName foo --variant:VCF " + testDir + "twoDeletions.vcf", 1,
                Arrays.asList("637256ee5348c1c57f1dadf581b06ed9"));
        executeTest("test deletions", spec);
    }
}
