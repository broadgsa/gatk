package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.util.Arrays;

public class VariantAnnotatorIntegrationTest extends WalkerTest {

    public static String baseTestString() {
        return "-T VariantAnnotator -R " + b36KGReference + " --no_cmdline_in_header -o %s";
    }

    @Test
    public void testHasAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + privateTestDir + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("360610e4990860bb5c45249b8ac31e5b"));
        executeTest("test file has annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsNotAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + privateTestDir + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("d69a3c92a0e8f44e09e7377e3eaed4e8"));
        executeTest("test file has annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testHasAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant " + privateTestDir + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("e0a08416249515ea18bd0663c90c9330"));
        executeTest("test file has annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant " + privateTestDir + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("0b60da46ba0eabb3abe5e0288937f9b0"));
        executeTest("test file has annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + privateTestDir + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("540a9be8a8cb85b0f675fea1184bf78c"));
        executeTest("test file doesn't have annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsNotAsking2() {
        // the genotype annotations in this file are actually out of order.  If you don't parse the genotypes
        // they don't get reordered.  It's a good test of the genotype ordering system.
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + privateTestDir + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("f900e65b65ff0f9d9eb0891ef9b28c73"));
        executeTest("test file doesn't have annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant " + privateTestDir + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("5eb576d0234c912d8efea184492691d0"));
        executeTest("test file doesn't have annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant " + privateTestDir + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("8860524d793d24b2e32f318433fcf527"));
        executeTest("test file doesn't have annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testExcludeAnnotations() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard -XA FisherStrand -XA ReadPosRankSumTest --variant " + privateTestDir + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("f33f417fad98c05d9cd08ffa22943b0f"));
        executeTest("test exclude annotations", spec);
    }

    @Test
    public void testOverwritingHeader() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant " + privateTestDir + "vcfexample4.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,001,292", 1,
                Arrays.asList("ebbf32f5b8b8d22f2eb247a0a3db3da0"));
        executeTest("test overwriting header", spec);
    }

    @Test
    public void testNoReads() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant " + privateTestDir + "vcfexample3empty.vcf -L " + privateTestDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("1c423b7730b9805e7b885ece924286e0"));
        executeTest("not passing it any reads", spec);
    }

    @Test
    public void testDBTagWithDbsnp() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --dbsnp " + b36dbSNP129 + " -G Standard --variant " + privateTestDir + "vcfexample3empty.vcf -L " + privateTestDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("54d7d5bb9404652857adf5e50d995f30"));
        executeTest("getting DB tag with dbSNP", spec);
    }

    @Test
    public void testMultipleIdsWithDbsnp() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --alwaysAppendDbsnpId --dbsnp " + b36dbSNP129 + " -G Standard --variant " + privateTestDir + "vcfexample3withIDs.vcf -L " + privateTestDir + "vcfexample3withIDs.vcf", 1,
                Arrays.asList("5fe63e511061ed4f91d938e72e7e3c39"));
        executeTest("adding multiple IDs with dbSNP", spec);
    }

    @Test
    public void testDBTagWithHapMap() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --comp:H3 " + privateTestDir + "fakeHM3.vcf -G Standard --variant " + privateTestDir + "vcfexample3empty.vcf -L " + privateTestDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("cc7184263975595a6e2473d153227146"));
        executeTest("getting DB tag with HM3", spec);
    }

    @Test
    public void testNoQuals() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + privateTestDir + "noQual.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L " + privateTestDir + "noQual.vcf -A QualByDepth", 1,
                Arrays.asList("aea983adc01cd059193538cc30adc17d"));
        executeTest("test file doesn't have QUALs", spec);
    }

    @Test
    public void testUsingExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --resource:foo " + privateTestDir + "targetAnnotations.vcf -G Standard --variant " + privateTestDir + "vcfexample3empty.vcf -E foo.AF -L " + privateTestDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("2b0e8cdfd691779befc5ac123d1a1887"));
        executeTest("using expression", spec);
    }

    @Test
    public void testUsingExpressionWithID() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --resource:foo " + privateTestDir + "targetAnnotations.vcf -G Standard --variant " + privateTestDir + "vcfexample3empty.vcf -E foo.ID -L " + privateTestDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("3de1d1998203518098ffae233f3e2352"));
        executeTest("using expression with ID", spec);
    }

    @Test
    public void testTabixAnnotations() {
        final String MD5 = "99938d1e197b8f10c408cac490a00a62";
        for ( String file : Arrays.asList("CEU.exon.2010_03.sites.vcf", "CEU.exon.2010_03.sites.vcf.gz")) {
            WalkerTestSpec spec = new WalkerTestSpec(
                    baseTestString() + " -A HomopolymerRun --variant:vcf " + validationDataLocation + file + " -L " + validationDataLocation + "CEU.exon.2010_03.sites.vcf --no_cmdline_in_header", 1,
                    Arrays.asList(MD5));
            executeTest("Testing lookup vcf tabix vs. vcf tribble", spec);
        }
    }

    @Test
    public void testSnpEffAnnotations() {
        WalkerTestSpec spec = new WalkerTestSpec(
            "-T VariantAnnotator -R " + hg19Reference + " --no_cmdline_in_header -o %s -A SnpEff --variant " +
            validationDataLocation + "1kg_exomes_unfiltered.AFR.unfiltered.vcf --snpEffFile  " + validationDataLocation +
            "snpEff2.0.5.AFR.unfiltered.vcf -L 1:1-1,500,000 -L 2:232,325,429",
            1,
            Arrays.asList("d9291845ce5a8576898d293a829a05b7")
        );
        executeTest("Testing SnpEff annotations", spec);
    }

    @Test
    public void testSnpEffAnnotationsUnsupportedVersion() {
        WalkerTestSpec spec = new WalkerTestSpec(
            "-T VariantAnnotator -R " + hg19Reference + " --no_cmdline_in_header -o %s -A SnpEff --variant " +
            validationDataLocation + "1kg_exomes_unfiltered.AFR.unfiltered.vcf --snpEffFile  " + validationDataLocation +
            "snpEff.AFR.unfiltered.unsupported.version.vcf -L 1:1-1,500,000",
            1,
            UserException.class
        );
        executeTest("Testing SnpEff annotations (unsupported version)", spec);
    }

    @Test
    public void testTDTAnnotation() {
        final String MD5 = "427dfdc665359b67eff210f909ebf8a2";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + b37KGReference + " -A TransmissionDisequilibriumTest --variant:vcf " + privateTestDir + "ug.random50000.subset300bp.chr1.family.vcf" +
                        " -L " + privateTestDir + "ug.random50000.subset300bp.chr1.family.vcf --no_cmdline_in_header -ped " + privateTestDir + "ug.random50000.family.ped -o %s", 1,
                Arrays.asList(MD5));
        executeTest("Testing TDT annotation ", spec);
    }


    @Test
    public void testChromosomeCountsPed() {
        final String MD5 = "6b5cbedf4a8b3385edf128d81c8a46f2";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + b37KGReference + " -A ChromosomeCounts --variant:vcf " + privateTestDir + "ug.random50000.subset300bp.chr1.family.vcf" +
                        " -L " + privateTestDir + "ug.random50000.subset300bp.chr1.family.vcf --no_cmdline_in_header -ped " + privateTestDir + "ug.random50000.family.ped -o %s", 1,
                Arrays.asList(MD5));
        executeTest("Testing ChromosomeCounts annotation with PED file", spec);
    }

    @Test
    public void testInbreedingCoeffPed() {
        final String MD5 = "159a771c1deaeffb786097e106943893";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + b37KGReference + " -A InbreedingCoeff --variant:vcf " + privateTestDir + "ug.random50000.subset300bp.chr1.family.vcf" +
                        " -L " + privateTestDir + "ug.random50000.subset300bp.chr1.family.vcf --no_cmdline_in_header -ped " + privateTestDir + "ug.random50000.family.ped -o %s", 1,
                Arrays.asList(MD5));
        executeTest("Testing InbreedingCoeff annotation with PED file", spec);
    }

}
