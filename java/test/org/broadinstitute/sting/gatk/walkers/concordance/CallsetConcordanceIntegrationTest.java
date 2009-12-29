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
                Arrays.asList("f4c4b7430f0e293c27ce38cb89b9338b"));
        executeTest("testSimpleVenn", spec);
    }

    @Test
    public void testSNPConcordance() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "NA12878.example1.vcf -B set2,VCF," + validationDataLocation + "NA12878.example2.vcf -CT SNPGenotypeConcordance:qscore=5", 1,
                Arrays.asList("f754e046bc0fa3a4b3430061e412ef0d"));
        executeTest("testSNPConcordance", spec);
    }

    @Test
    public void testNWayVenn() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "NA12878.example1.vcf -B set2,VCF," + validationDataLocation + "NA12878.example2.vcf -B set3,VCF," + validationDataLocation + "CEU.sample.vcf -CT NWayVenn", 1,
                Arrays.asList("0527ea8ec7de3a144bd0a56db80d62ba"));
        executeTest("testNWayVenn", spec);
    }

    @Test
    public void testMulti() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF," + validationDataLocation + "NA12878.example1.vcf -B set2,VCF," + validationDataLocation + "NA12878.example2.vcf -CT SimpleVenn -CT NWayVenn -CT SNPGenotypeConcordance:qscore=5", 1,
                Arrays.asList("6fbe00cb68d2cdc59dfcb79024fd9893"));
        executeTest("testMulti", spec);
    }
}