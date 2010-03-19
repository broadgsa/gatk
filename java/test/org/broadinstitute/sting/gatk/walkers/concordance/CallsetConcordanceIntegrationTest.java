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
                Arrays.asList("d9124d2b0fb5bec5bc50c26a16b4e900"));
        executeTest("testSimpleVenn", spec);
    }

    @Test
    public void testSNPConcordance() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "NA12878.example1.vcf -B set2,VCF," + validationDataLocation + "NA12878.example2.vcf -CT SNPGenotypeConcordance:qscore=5", 1,
                Arrays.asList("df1fbc744947f316f65f51a21368b0e4"));
        executeTest("testSNPConcordance", spec);
    }

    @Test
    public void testNWayVenn() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "NA12878.example1.vcf -B set2,VCF," + validationDataLocation + "NA12878.example2.vcf -B set3,VCF," + validationDataLocation + "CEU.sample.vcf -CT NWayVenn", 1,
                Arrays.asList("aa835ae5368b35f376b844d7f8ef2976"));
        executeTest("testNWayVenn", spec);
    }

    @Test
    public void testMulti() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "NA12878.example1.vcf -B set2,VCF," + validationDataLocation + "NA12878.example2.vcf -CT SimpleVenn -CT NWayVenn -CT SNPGenotypeConcordance:qscore=5", 1,
                Arrays.asList("c9ef68cc3b7dc08f1d2b49170e6560ab"));
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
