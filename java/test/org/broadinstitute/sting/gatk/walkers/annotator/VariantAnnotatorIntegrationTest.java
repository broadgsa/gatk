package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.StingException;
import org.junit.Test;

import java.util.Arrays;

public class VariantAnnotatorIntegrationTest extends WalkerTest {

    public static String secondBaseTestString() {
        return "-T VariantAnnotator -R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -vcf %s -A SecondBaseSkew";
    }

    public static String validationDataPath() {
         return validationDataLocation + "";
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
        return "-T VariantAnnotator -R " + oneKGLocation + "reference/human_b36_both.fasta -o %s";
    }



    @Test
    public void testHasAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("5a91f6b50fc136d7d3f8735dbc64defe"));
        executeTest("test file has annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsNotAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("0a561b161a06e68b88417ff5fe365871"));
        executeTest("test file has annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testHasAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G \"Standard\" -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("dff7c888fbd142b6eff27c3233c1292d"));
        executeTest("test file has annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G \"Standard\" -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("45c9d5136900d9ab347f489d6f442bb4"));
        executeTest("test file has annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("16e578597eed68609268c00886f2842a"));
        executeTest("test file doesn't have annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsNotAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("0e947087c44479faefe616f6b6f7d272"));
        executeTest("test file doesn't have annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G \"Standard\" -B variant,VCF," + validationDataLocation + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("7c13af3226a92bfc3826c37b30f179e0"));
        executeTest("test file doesn't have annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G \"Standard\" -B variant,VCF," + validationDataLocation + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("369ad538eacdc0943b58a807d92d823c"));
        executeTest("test file doesn't have annotations, asking for annotations, #2", spec);
    }

}