package org.broadinstitute.sting.gatk.walkers.concordance;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class CallsetConcordanceIntegrationTest extends WalkerTest {
    public static String baseTestString() {
        return "-T CallsetConcordance -R " + oneKGLocation + "reference/human_b36_both.fasta -L 1:1-8000 -CO %s";
    }

    @Test
    public void testSimpleVenn() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "NA12878.example1.vcf -B set2,VCF," + validationDataLocation + "NA12878.example2.vcf -CT SimpleVenn", 1,
                Arrays.asList("cdd027a0d7bbfae2ba75480fcaa14356"));
        executeTest("testSimpleVenn", spec);
    }

    @Test
    public void testSNPConcordance() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "NA12878.example1.vcf -B set2,VCF," + validationDataLocation + "NA12878.example2.vcf -CT SNPGenotypeConcordance:qscore=5", 1,
                Arrays.asList("310bf0824fc407b2e40cbbaea234471e"));
        executeTest("testSNPConcordance", spec);
    }

    @Test
    public void testNWayVenn() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "NA12878.example1.vcf -B set2,VCF," + validationDataLocation + "NA12878.example2.vcf -B set3,VCF," + validationDataLocation + "CEU.sample.vcf -CT NWayVenn", 1,
                Arrays.asList("179adafae6efdc879fc22442ab9e599f"));
        executeTest("testNWayVenn", spec);
    }

    @Test
    public void testMulti() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "NA12878.example1.vcf -B set2,VCF," + validationDataLocation + "NA12878.example2.vcf -CT SimpleVenn -CT NWayVenn -CT SNPGenotypeConcordance:qscore=5", 1,
                Arrays.asList("c8fe63633ef6ed2b6068958f06cddfe0"));
        executeTest("testMulti", spec);
    }

    @Test
    public void testComplex() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "complexExample.vcf -B set2,VCF," + validationDataLocation + "complexExample.vcf -CT NWayVenn", 1,
                Arrays.asList("7bc72ec5f8b0fda5d59ebd2526b53e48"));
        executeTest("testComplex", spec);
    }
}