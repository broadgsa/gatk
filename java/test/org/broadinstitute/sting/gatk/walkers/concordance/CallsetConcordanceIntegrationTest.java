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
                Arrays.asList("2c7e18901dbf27bac9f36b3dbee063c6"));
        executeTest("testSimpleVenn", spec);
    }

    @Test
    public void testSNPConcordance() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example1.vcf -B set2,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example2.vcf -CT SNPGenotypeConcordance:qscore=5", 1,
                Arrays.asList("c21d59fc3194c39c662d2e74b53dcf9c"));
        executeTest("testSNPConcordance", spec);
    }

    @Test
    public void testNWayVenn() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example1.vcf -B set2,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example2.vcf -B set3,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/CEU.sample.vcf -CT NWayVenn", 1,
                Arrays.asList("2b38ae235edd10773dbee0bfae036e35"));
        executeTest("testNWayVenn", spec);
    }

    @Test
    public void testMulti() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example1.vcf -B set2,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example2.vcf -CT SimpleVenn -CT NWayVenn -CT SNPGenotypeConcordance:qscore=5", 1,
                Arrays.asList("9bcc83aadac00a160cef20e7126368ee"));
        executeTest("testMulti", spec);
    }
}