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
                Arrays.asList("851b68004874f3a2e76d795e7401f8a0"));
        executeTest("testSimpleVenn", spec);
    }

    @Test
    public void testSNPConcordance() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example1.vcf -B set2,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example2.vcf -CT SNPGenotypeConcordance:qscore=5", 1,
                Arrays.asList("7afb56b30257fe2d66bee7a029d75685"));
        executeTest("testSNPConcordance", spec);
    }

    @Test
    public void testNWayVenn() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example1.vcf -B set2,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example2.vcf -B set3,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/CEU.sample.vcf -CT NWayVenn", 1,
                Arrays.asList("f452c04c600ad10c054f18b0c77b53d5"));
        executeTest("testNWayVenn", spec);
    }
}