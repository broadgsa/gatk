package org.broadinstitute.sting.gatk.walkers.concordance;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class CallsetConcordanceIntegrationTest extends WalkerTest {
    public static String baseTestString() {
        return "-T CallsetConcordance -R /broad/1KG/reference/human_b36_both.fasta -L 1:1-8000 -CO %s";
    }

    @Test
    public void testSimpleVenn() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example1.vcf -B set2,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example2.vcf -CT SimpleVenn", 1,
                Arrays.asList("0a71c8f06b4179ba59cefad962cd034c"));
        executeTest("testSimpleVenn", spec);
    }

    @Test
    public void testSNPConcordance() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example1.vcf -B set2,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example2.vcf -CT SNPGenotypeConcordance:qscore=5", 1,
                Arrays.asList("5da8bf664813f0ab8b22070097f6900e"));
        executeTest("testSNPConcordance", spec);
    }

    @Test
    public void testNWayVenn() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example1.vcf -B set2,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example2.vcf -B set3,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/CEU.sample.vcf -CT NWayVenn", 1,
                Arrays.asList("9da88442eea094da8b6110d8f5ed4408"));
        executeTest("testNWayVenn", spec);
    }
}