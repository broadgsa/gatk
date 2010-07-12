package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class VariantAnnotatorIntegrationTest extends WalkerTest {

    public static String baseTestString() {
        return "-T VariantAnnotator -R " + oneKGLocation + "reference/human_b36_both.fasta -o %s";
    }

    @Test
    public void testHasAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("92369e62b3bade625992bc1741e9367a"));
        executeTest("test file has annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsNotAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("1a394cdace4ff7698b603c2252ee9365"));
        executeTest("test file has annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testHasAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G \"Standard\" -B variant,VCF," + validationDataLocation + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("5a8c9c18da89b052e7f3797cd391bff6"));
        executeTest("test file has annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G \"Standard\" -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("0c8607c69d761bd058b944d1e98708b6"));
        executeTest("test file has annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("006642aabd41e806d3ea64bf5441f1d2"));
        executeTest("test file doesn't have annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsNotAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("c66ed83d2a981c43546d05120809ed57"));
        executeTest("test file doesn't have annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G \"Standard\" -B variant,VCF," + validationDataLocation + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("3b570e13bd0e1cb5750da268178ce656"));
        executeTest("test file doesn't have annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G \"Standard\" -B variant,VCF," + validationDataLocation + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("42e9c57cd5d266fd7837354fc8495176"));
        executeTest("test file doesn't have annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testNoReads() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G \"Standard\" -B variant,VCF," + validationDataLocation + "vcfexample3empty.vcf -BTI variant", 1,
                Arrays.asList("374da32e1b9540bfe8ee6ea5b896babd"));
        executeTest("not passing it any reads", spec);
    }

    @Test
    public void testDBTag() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -D " + GATKDataLocation + "dbsnp_129_b36.rod -G \"Standard\" -B variant,VCF," + validationDataLocation + "vcfexample3empty.vcf -BTI variant", 1,
                Arrays.asList("39a6b268e3edc1eee5f71d77460e448e"));
        executeTest("getting DB tag", spec);
    }
}
