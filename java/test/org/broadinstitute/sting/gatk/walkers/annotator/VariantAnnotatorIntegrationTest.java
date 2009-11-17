package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class VariantAnnotatorIntegrationTest extends WalkerTest {
    public static String baseTestString() {
        return "-T VariantAnnotator -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,010,000 -vcf %s";
    }

    @Test
    public void testNoAnnots1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf ", 1,
                Arrays.asList("4e231f16c202d88ca3adb17168af0e0f"));
        //executeTest("testNoAnnots1", spec);
    }

    @Test
    public void testNoAnnots2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample3.vcf ", 1,
                Arrays.asList("ef0c59e47a2afcbecf2bcef6aa01e7e7"));
        //executeTest("testNoAnnots2", spec);
    }

    @Test
    public void testAllAnnots1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -all -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf ", 1,
                Arrays.asList("ced92b5ac9e2692c4d8acce1235317b6"));
        //executeTest("testAllAnnots1", spec);
    }

    @Test
    public void testAllAnnots2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -all -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample3.vcf ", 1,
                Arrays.asList("573a6c02f659b2c4cf014f84bd0b9c8a"));
        //executeTest("testAllAnnots2", spec);
    }
}