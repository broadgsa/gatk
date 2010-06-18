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
                Arrays.asList("a1970effe9c51923d52af9034e778de4"));
        executeTest("testSimpleVenn", spec);
    }

    @Test
    public void testSNPConcordance() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "NA12878.example1.vcf -B set2,VCF," + validationDataLocation + "NA12878.example2.vcf -CT SNPGenotypeConcordance:qscore=5", 1,
                Arrays.asList("e7a0d52c266ba3c76283111674c7168f"));
        executeTest("testSNPConcordance", spec);
    }

    @Test
    public void testNWayVenn() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "NA12878.example1.vcf -B set2,VCF," + validationDataLocation + "NA12878.example2.vcf -B set3,VCF," + validationDataLocation + "CEU.sample.vcf -CT NWayVenn", 1,
                Arrays.asList("e65fc811137fca7d6c32125240c7468f"));
        executeTest("testNWayVenn", spec);
    }

    @Test
    public void testMulti() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "NA12878.example1.vcf -B set2,VCF," + validationDataLocation + "NA12878.example2.vcf -CT SimpleVenn -CT NWayVenn -CT SNPGenotypeConcordance:qscore=5", 1,
                Arrays.asList("ddc2507590e28743e9cb4b132cb066e7"));
        executeTest("testMulti", spec);
    }

    @Test
    public void testComplex() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "complexExample.vcf -B set2,VCF," + validationDataLocation + "complexExample.vcf -CT NWayVenn", 1,
                Arrays.asList("250df7bde7a8cf9c7ee7c5704183ea88"));
        executeTest("testComplex", spec);
    }
}
