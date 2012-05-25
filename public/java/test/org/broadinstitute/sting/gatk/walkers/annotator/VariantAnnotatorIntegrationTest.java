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
                baseTestString() + " --variant:VCF3 " + testDir + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("8a105fa5eebdfffe7326bc5b3d8ffd1c"));
        executeTest("test file has annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsNotAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant:VCF3 " + testDir + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("964f1016ec9a3c55333f62dd834c14d6"));
        executeTest("test file has annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testHasAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant:VCF3 " + testDir + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("21810e5b51a0431b510a8d7a73f8771b"));
        executeTest("test file has annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant:VCF3 " + testDir + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("31263a127c606f72fb00fac69f0320a6"));
        executeTest("test file has annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant:VCF3 " + testDir + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("42ccee09fa9f8c58f4a0d4f1139c094f"));
        executeTest("test file doesn't have annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsNotAsking2() {
        // the genotype annotations in this file are actually out of order.  If you don't parse the genotypes
        // they don't get reordered.  It's a good test of the genotype ordering system.
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant:VCF3 " + testDir + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("f2ddfa8105c290b1f34b7a261a02a1ac"));
        executeTest("test file doesn't have annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant:VCF3 " + testDir + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("bd367dbcbcc45cc43588e1d7309a7ad1"));
        executeTest("test file doesn't have annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant:VCF3 " + testDir + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("7f200756227a1ee75f33c466955015c5"));
        executeTest("test file doesn't have annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testExcludeAnnotations() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard -XA FisherStrand -XA ReadPosRankSumTest --variant:VCF3 " + testDir + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("65dc85869239d52e3893e5ba47c1afac"));
        executeTest("test exclude annotations", spec);
    }

    @Test
    public void testOverwritingHeader() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant " + testDir + "vcfexample4.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,001,292", 1,
                Arrays.asList("1d98be77dad9c703402de0315db5176a"));
        executeTest("test overwriting header", spec);
    }

    @Test
    public void testNoReads() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant " + testDir + "vcfexample3empty.vcf -L " + testDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("920d796cfaa19c20153181d7be627435"));
        executeTest("not passing it any reads", spec);
    }

    @Test
    public void testDBTagWithDbsnp() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --dbsnp " + b36dbSNP129 + " -G Standard --variant " + testDir + "vcfexample3empty.vcf -L " + testDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("3a1e3ddd8c563a9de522c551e3d7e724"));
        executeTest("getting DB tag with dbSNP", spec);
    }

    @Test
    public void testMultipleIdsWithDbsnp() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --alwaysAppendDbsnpId --dbsnp " + b36dbSNP129 + " -G Standard --variant " + testDir + "vcfexample3withIDs.vcf -L " + testDir + "vcfexample3withIDs.vcf", 1,
                Arrays.asList("06107c9dfefc1ab53e366d0f86502ff2"));
        executeTest("adding multiple IDs with dbSNP", spec);
    }

    @Test
    public void testDBTagWithHapMap() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --comp:H3 " + testDir + "fakeHM3.vcf -G Standard --variant " + testDir + "vcfexample3empty.vcf -L " + testDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("fe1e76a7a0d95ee7c604194b165a640d"));
        executeTest("getting DB tag with HM3", spec);
    }

    @Test
    public void testNoQuals() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + testDir + "noQual.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L " + testDir + "noQual.vcf -A QualByDepth", 1,
                Arrays.asList("e531c9f90c17f0f859cd1ac851a8edd8"));
        executeTest("test file doesn't have QUALs", spec);
    }

    @Test
    public void testUsingExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --resource:foo " + testDir + "targetAnnotations.vcf -G Standard --variant:VCF3 " + testDir + "vcfexample3empty.vcf -E foo.AF -L " + testDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("b5583188347922476e23080a265c214e"));
        executeTest("using expression", spec);
    }

    @Test
    public void testUsingExpressionWithID() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --resource:foo " + testDir + "targetAnnotations.vcf -G Standard --variant:VCF3 " + testDir + "vcfexample3empty.vcf -E foo.ID -L " + testDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("bb6d9553d6954923cdef157753735b87"));
        executeTest("using expression with ID", spec);
    }

    @Test
    public void testTabixAnnotations() {
        final String MD5 = "bb9a148716fc69d706c5be146c1afa00";
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
            Arrays.asList("ffbda45b3682c9b83cb541d83f6c15d6")
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
        final String MD5 = "a78c1e950740d3c13c0258960c5fa8e1";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + b37KGReference + " -A TransmissionDisequilibriumTest --variant:vcf " + testDir + "ug.random50000.subset300bp.chr1.family.vcf" +
                        " -L " + testDir + "ug.random50000.subset300bp.chr1.family.vcf --no_cmdline_in_header -ped " + testDir + "ug.random50000.family.ped -o %s", 1,
                Arrays.asList(MD5));
        executeTest("Testing TDT annotation ", spec);
    }


    @Test
    public void testChromosomeCountsPed() {
        final String MD5 = "32df3ceb63c277df442ed55fb8684933";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + b37KGReference + " -A ChromosomeCounts --variant:vcf " + testDir + "ug.random50000.subset300bp.chr1.family.vcf" +
                        " -L " + testDir + "ug.random50000.subset300bp.chr1.family.vcf --no_cmdline_in_header -ped " + testDir + "ug.random50000.family.ped -o %s", 1,
                Arrays.asList(MD5));
        executeTest("Testing ChromosomeCounts annotation with PED file", spec);
    }

    @Test
    public void testInbreedingCoeffPed() {
        final String MD5 = "7f1314fada5cb1f35ba1996f8a7a686b";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + b37KGReference + " -A InbreedingCoeff --variant:vcf " + testDir + "ug.random50000.subset300bp.chr1.family.vcf" +
                        " -L " + testDir + "ug.random50000.subset300bp.chr1.family.vcf --no_cmdline_in_header -ped " + testDir + "ug.random50000.family.ped -o %s", 1,
                Arrays.asList(MD5));
        executeTest("Testing InbreedingCoeff annotation with PED file", spec);
    }

}
