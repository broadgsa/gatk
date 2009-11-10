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
                Arrays.asList("1b8e26cd30e993da9318abd6475f38d0"));
        executeTest("testSimpleVenn", spec);
    }

    @Test
    public void testSNPConcordance() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example1.vcf -B set2,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example2.vcf -CT SNPGenotypeConcordance:qscore=5", 1,
                Arrays.asList("5a89b8edcdf2e3f469ac354cb1524033"));
        executeTest("testSNPConcordance", spec);
    }

    @Test
    public void testNWayVenn() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B set1,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example1.vcf -B set2,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.example2.vcf -B set3,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/CEU.sample.vcf -CT NWayVenn", 1,
                Arrays.asList("1dec083580b75a9c59fcb61426117134"));
        executeTest("testNWayVenn", spec);
    }
}