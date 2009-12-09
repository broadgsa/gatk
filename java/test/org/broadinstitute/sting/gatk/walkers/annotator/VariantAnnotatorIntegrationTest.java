package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.StingException;
import org.junit.Test;

import java.util.Arrays;

public class VariantAnnotatorIntegrationTest extends WalkerTest {

    public static String secondBaseTestString() {
        return "-T VariantAnnotator -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -vcf %s -A SecondBaseSkew";
    }

    public static String validationDataPath() {
         return "/humgen/gsa-scr1/GATK_Data/Validation_Data/";
    }

    public static String secondBaseTestFile( int testNo ) {
        switch ( testNo ) {
            case 1: return "NA12891";
            case 2: return "NA20762";
            default: throw new StingException("Impossible test has been run: secondbasetest number "+ testNo);
        }
    }

    public static String secondBaseTestInterval ( int testNo ) {
        switch ( testNo ) {
            case 1: return "-L chr1:14,000,000-18,000,000";
            case 2: return "-L chr22:20660081-20660083 -L chr22:29198100-29198104 -L chr22:29821330-29821334";
            default: throw new StingException("Impossible test has been run: secondbasetest number "+testNo);
        }
    }

    public static String secondBaseTestmd5( int testNo ) {
        switch ( testNo ) {
            case 1: return "bf64bac186fd682018dd7f0419d90190";
            case 2: return "67f40627b12be31efe02c9d853fbcf37";
            default: throw new StingException("Impossible test has been run: secondbasetest number "+testNo);
        }
    }
    
    public static String baseTestString() {
        return "-T VariantAnnotator -R /broad/1KG/reference/human_b36_both.fasta -vcf %s";
    }



    @Test
    public void testHasAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf -I /humgen/gsa-scr1/GATK_Data/Validation_Data/low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("6903f8b31820ce7a56230b99f9a9309c"));
        executeTest("test file has annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsNotAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample3.vcf -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("5298c24e956361d209f14ac6138a3bbd"));
        executeTest("test file has annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testHasAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -standard -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2.vcf -I /humgen/gsa-scr1/GATK_Data/Validation_Data/low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("24609ede968a90ebff949791694f70a0"));
        executeTest("test file has annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -standard -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample3.vcf -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("9f05c78ef9913259c02f8bda6ec97505"));
        executeTest("test file has annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2empty.vcf -I /humgen/gsa-scr1/GATK_Data/Validation_Data/low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("b46f395864cb71b887d69342f56e7fdb"));
        executeTest("test file doesn't have annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsNotAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample3empty.vcf -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("aaae89c48fcab615fe4204220ec62859"));
        executeTest("test file doesn't have annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -standard -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample2empty.vcf -I /humgen/gsa-scr1/GATK_Data/Validation_Data/low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("54811556d642a83fcfa2b359aa275023"));
        executeTest("test file doesn't have annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -standard -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample3empty.vcf -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("b67282a03c92ae3d24a7b588aa1e61df"));
        executeTest("test file doesn't have annotations, asking for annotations, #2", spec);
    }

}