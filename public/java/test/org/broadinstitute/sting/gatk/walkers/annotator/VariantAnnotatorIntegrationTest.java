package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.util.Arrays;

public class VariantAnnotatorIntegrationTest extends WalkerTest {

    public static String baseTestString() {
        return "-T VariantAnnotator -R " + b36KGReference + " -NO_HEADER -o %s";
    }

    @Test
    public void testHasAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("8a105fa5eebdfffe7326bc5b3d8ffd1c"));
        executeTest("test file has annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsNotAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant:VCF3 " + validationDataLocation + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("964f1016ec9a3c55333f62dd834c14d6"));
        executeTest("test file has annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testHasAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("8e7de435105499cd71ffc099e268a83e"));
        executeTest("test file has annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant:VCF3 " + validationDataLocation + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("64b6804cb1e27826e3a47089349be581"));
        executeTest("test file has annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant:VCF3 " + validationDataLocation + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("42ccee09fa9f8c58f4a0d4f1139c094f"));
        executeTest("test file doesn't have annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsNotAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant:VCF3 " + validationDataLocation + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("f2ddfa8105c290b1f34b7a261a02a1ac"));
        executeTest("test file doesn't have annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant:VCF3 " + validationDataLocation + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("fd1ffb669800c2e07df1e2719aa38e49"));
        executeTest("test file doesn't have annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant:VCF3 " + validationDataLocation + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("09f8e840770a9411ff77508e0ed0837f"));
        executeTest("test file doesn't have annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testExcludeAnnotations() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard -XA FisherStrand -XA ReadPosRankSumTest --variant:VCF3 " + validationDataLocation + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("b49fe03aa4b675db80a9db38a3552c95"));
        executeTest("test exclude annotations", spec);
    }

    @Test
    public void testOverwritingHeader() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant " + validationDataLocation + "vcfexample4.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,001,292", 1,
                Arrays.asList("78d2c19f8107d865970dbaf3e12edd92"));
        executeTest("test overwriting header", spec);
    }

    @Test
    public void testNoReads() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant " + validationDataLocation + "vcfexample3empty.vcf -L " + validationDataLocation + "vcfexample3empty.vcf", 1,
                Arrays.asList("16e3a1403fc376320d7c69492cad9345"));
        executeTest("not passing it any reads", spec);
    }

    @Test
    public void testDBTagWithDbsnp() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --dbsnp " + b36dbSNP129 + " -G Standard --variant " + validationDataLocation + "vcfexample3empty.vcf -L " + validationDataLocation + "vcfexample3empty.vcf", 1,
                Arrays.asList("3da8ca2b6bdaf6e92d94a8c77a71313d"));
        executeTest("getting DB tag with dbSNP", spec);
    }

    @Test
    public void testDBTagWithHapMap() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --comp:H3 " + validationDataLocation + "fakeHM3.vcf -G Standard --variant " + validationDataLocation + "vcfexample3empty.vcf -L " + validationDataLocation + "vcfexample3empty.vcf", 1,
                Arrays.asList("1bc01c5b3bd0b7aef75230310c3ce688"));
        executeTest("getting DB tag with HM3", spec);
    }

    @Test
    public void testUsingExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --resource:foo " + validationDataLocation + "targetAnnotations.vcf -G Standard --variant:VCF3 " + validationDataLocation + "vcfexample3empty.vcf -E foo.AF -L " + validationDataLocation + "vcfexample3empty.vcf", 1,
                Arrays.asList("ae30a1ac7bfbc3d22a327f8b689cad31"));
        executeTest("using expression", spec);
    }

    @Test
    public void testTabixAnnotations() {
        final String MD5 = "13269d5a2e16f06fd755cc0fb9271acf";
        for ( String file : Arrays.asList("CEU.exon.2010_03.sites.vcf", "CEU.exon.2010_03.sites.vcf.gz")) {
            WalkerTestSpec spec = new WalkerTestSpec(
                    baseTestString() + " -A HomopolymerRun --variant:vcf " + validationDataLocation + file + " -L " + validationDataLocation + "CEU.exon.2010_03.sites.vcf -NO_HEADER", 1,
                    Arrays.asList(MD5));
            executeTest("Testing lookup vcf tabix vs. vcf tribble", spec);
        }
    }

    @Test
    public void testSnpEffAnnotations() {
        WalkerTestSpec spec = new WalkerTestSpec(
            "-T VariantAnnotator -R " + hg19Reference + " -NO_HEADER -o %s -A SnpEff --variant " +
            validationDataLocation + "1kg_exomes_unfiltered.AFR.unfiltered.vcf --snpEffFile  " + validationDataLocation +
            "snpEff.AFR.unfiltered.vcf -L 1:1-1,500,000 -L 2:232,325,429",
            1,
            Arrays.asList("122321a85e448f21679f6ca15c5e22ad")
        );
        executeTest("Testing SnpEff annotations", spec);
    }

    @Test
    public void testSnpEffAnnotationsUnsupportedVersion() {
        WalkerTestSpec spec = new WalkerTestSpec(
            "-T VariantAnnotator -R " + hg19Reference + " -NO_HEADER -o %s -A SnpEff --variant " +
            validationDataLocation + "1kg_exomes_unfiltered.AFR.unfiltered.vcf --snpEffFile  " + validationDataLocation +
            "snpEff.AFR.unfiltered.unsupported.version.vcf -L 1:1-1,500,000",
            1,
            UserException.class
        );
        executeTest("Testing SnpEff annotations (unsupported version)", spec);
    }
}
