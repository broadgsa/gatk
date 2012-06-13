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
                Arrays.asList("dfa5dff09fa964b06da19c0f4aff6928"));
        executeTest("test file has annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsNotAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant:VCF3 " + testDir + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("9914bd19f6235c550e5182e0f4591da6"));
        executeTest("test file has annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testHasAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant:VCF3 " + testDir + "vcfexample2.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("6a52ef10bb10d72cdd82a8f7afc2dd09"));
        executeTest("test file has annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testHasAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant:VCF3 " + testDir + "vcfexample3.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("edf35b2cfd033d5abafc33fe739759e2"));
        executeTest("test file has annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsNotAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant:VCF3 " + testDir + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("dd89dfa22f0e1d6760095e04f528d62a"));
        executeTest("test file doesn't have annotations, not asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsNotAsking2() {
        // the genotype annotations in this file are actually out of order.  If you don't parse the genotypes
        // they don't get reordered.  It's a good test of the genotype ordering system.
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant:VCF3 " + testDir + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("542d9ed8290ef7868387af4127e0b5fa"));
        executeTest("test file doesn't have annotations, not asking for annotations, #2", spec);
    }

    @Test
    public void testNoAnnotsAsking1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant:VCF3 " + testDir + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("b1b32ed3b831c92c94258c8e4a60e8c9"));
        executeTest("test file doesn't have annotations, asking for annotations, #1", spec);
    }

    @Test
    public void testNoAnnotsAsking2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant:VCF3 " + testDir + "vcfexample3empty.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000", 1,
                Arrays.asList("7e5b8ad0c7005036e39739de77546023"));
        executeTest("test file doesn't have annotations, asking for annotations, #2", spec);
    }

    @Test
    public void testExcludeAnnotations() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard -XA FisherStrand -XA ReadPosRankSumTest --variant:VCF3 " + testDir + "vcfexample2empty.vcf -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("ef046909a6f6c6cb43653a255a99a014"));
        executeTest("test exclude annotations", spec);
    }

    @Test
    public void testOverwritingHeader() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant " + testDir + "vcfexample4.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,001,292", 1,
                Arrays.asList("341fbac1802c54b1ee658b4bb54e963e"));
        executeTest("test overwriting header", spec);
    }

    @Test
    public void testNoReads() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -G Standard --variant " + testDir + "vcfexample3empty.vcf -L " + testDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("f44f1a10a2e3b8e749d1903db7429f44"));
        executeTest("not passing it any reads", spec);
    }

    @Test
    public void testDBTagWithDbsnp() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --dbsnp " + b36dbSNP129 + " -G Standard --variant " + testDir + "vcfexample3empty.vcf -L " + testDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("cb0ee4261b72caf652c77c6d75bc7c85"));
        executeTest("getting DB tag with dbSNP", spec);
    }

    @Test
    public void testMultipleIdsWithDbsnp() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --alwaysAppendDbsnpId --dbsnp " + b36dbSNP129 + " -G Standard --variant " + testDir + "vcfexample3withIDs.vcf -L " + testDir + "vcfexample3withIDs.vcf", 1,
                Arrays.asList("e2b9cb00137bb4f2cb86b77d30b35dd3"));
        executeTest("adding multiple IDs with dbSNP", spec);
    }

    @Test
    public void testDBTagWithHapMap() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --comp:H3 " + testDir + "fakeHM3.vcf -G Standard --variant " + testDir + "vcfexample3empty.vcf -L " + testDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("b501a8f92abb1623e71c4593b68e8b31"));
        executeTest("getting DB tag with HM3", spec);
    }

    @Test
    public void testNoQuals() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant " + testDir + "noQual.vcf -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L " + testDir + "noQual.vcf -A QualByDepth", 1,
                Arrays.asList("4a247f039dfb16ac05b38a0dd5f98da6"));
        executeTest("test file doesn't have QUALs", spec);
    }

    @Test
    public void testUsingExpression() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --resource:foo " + testDir + "targetAnnotations.vcf -G Standard --variant:VCF3 " + testDir + "vcfexample3empty.vcf -E foo.AF -L " + testDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("33926902898f1c9ad871afe54c64e5b5"));
        executeTest("using expression", spec);
    }

    @Test
    public void testUsingExpressionWithID() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --resource:foo " + testDir + "targetAnnotations.vcf -G Standard --variant:VCF3 " + testDir + "vcfexample3empty.vcf -E foo.ID -L " + testDir + "vcfexample3empty.vcf", 1,
                Arrays.asList("636d9557a04648e639674069004a2315"));
        executeTest("using expression with ID", spec);
    }

    @Test
    public void testTabixAnnotations() {
        final String MD5 = "5aebcf8f76c649d645708b1262185c80";
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
            Arrays.asList("0c20cda1cf0b903a287f1807ae5bee02")
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
        final String MD5 = "240e31b14db09b5ce43c40e21ca52dc2";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + b37KGReference + " -A TransmissionDisequilibriumTest --variant:vcf " + testDir + "ug.random50000.subset300bp.chr1.family.vcf" +
                        " -L " + testDir + "ug.random50000.subset300bp.chr1.family.vcf --no_cmdline_in_header -ped " + testDir + "ug.random50000.family.ped -o %s", 1,
                Arrays.asList(MD5));
        executeTest("Testing TDT annotation ", spec);
    }


    @Test
    public void testChromosomeCountsPed() {
        final String MD5 = "d8ff51b00001e6f1b2cb7347b024c8bf";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + b37KGReference + " -A ChromosomeCounts --variant:vcf " + testDir + "ug.random50000.subset300bp.chr1.family.vcf" +
                        " -L " + testDir + "ug.random50000.subset300bp.chr1.family.vcf --no_cmdline_in_header -ped " + testDir + "ug.random50000.family.ped -o %s", 1,
                Arrays.asList(MD5));
        executeTest("Testing ChromosomeCounts annotation with PED file", spec);
    }

    @Test
    public void testInbreedingCoeffPed() {
        final String MD5 = "e94d589b5691e3ecfd9cc9475a384890";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T VariantAnnotator -R " + b37KGReference + " -A InbreedingCoeff --variant:vcf " + testDir + "ug.random50000.subset300bp.chr1.family.vcf" +
                        " -L " + testDir + "ug.random50000.subset300bp.chr1.family.vcf --no_cmdline_in_header -ped " + testDir + "ug.random50000.family.ped -o %s", 1,
                Arrays.asList(MD5));
        executeTest("Testing InbreedingCoeff annotation with PED file", spec);
    }

}
