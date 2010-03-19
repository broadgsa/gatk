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
                Arrays.asList("af14c7d03e3516b61f2702c9e4a7780f"));
        executeTest("test no action", spec);
    }

    @Test
    public void testClusteredSnps() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -window 10 -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("b2b18929289bb47f07bcc23d4cec94c4"));
        executeTest("test clustered SNPs", spec);
    }

    @Test
    public void testMask() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -mask foo -B mask,VCF," + validationDataLocation + "vcfexample2.vcf -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("af048053eaef84fc4d61c51c50be1e0a"));
        executeTest("test mask", spec);
    }

    @Test
    public void testFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("43a580d9c684a496f21d0d42939dd910"));
        executeTest("test filter #1", spec);
    }

    @Test
    public void testFilter2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("2bd1e8975d8105dec2fb6055fbf00569"));
        executeTest("test filter #2", spec);
    }

    @Test
    public void testFilterWithSeparateNames() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterName ABF -filter 'AlleleBalance < 70.0' --filterName FSF -filter 'FisherStrand == 1.4' -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("3479158ffd02a45371b5103277a30a53"));
        executeTest("test filter with separate names #2", spec);
    }
}