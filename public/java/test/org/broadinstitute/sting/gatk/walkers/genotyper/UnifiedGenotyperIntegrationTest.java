package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

// ********************************************************************************** //
// Note that this class also serves as an integration test for the VariantAnnotator!  //
// ********************************************************************************** //

public class UnifiedGenotyperIntegrationTest extends WalkerTest {

    private final static String baseCommand = "-T UnifiedGenotyper -R " + b36KGReference + " -nosl --no_cmdline_in_header -glm BOTH -minIndelFrac 0.0 --dbsnp " + b36dbSNP129;
    private final static String baseCommandIndels = "-T UnifiedGenotyper -R " + b36KGReference + " -nosl --no_cmdline_in_header -glm INDEL -mbq 20 -minIndelFrac 0.0 --dbsnp " + b36dbSNP129;
    private final static String baseCommandIndelsb37 = "-T UnifiedGenotyper -R " + b37KGReference + " -nosl --no_cmdline_in_header -glm INDEL -mbq 20 --dbsnp " + b37dbSNP132;

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing normal calling
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testMultiSamplePilot1() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10,022,000-10,025,000", 1,
                Arrays.asList("40dae29480e26efb117344e0e82e997c"));
        executeTest("test MultiSample Pilot1", spec);
    }

    @Test
    public void testWithAllelesPassedIn1() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + testDir + "allelesForUG.vcf -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,025,000", 1,
                Arrays.asList("96b9f6005d9437bf527344cf98319739"));
        executeTest("test MultiSample Pilot2 with alleles passed in", spec1);
    }

    @Test
    public void testWithAllelesPassedIn2() {
        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + testDir + "allelesForUG.vcf -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,025,000", 1,
                Arrays.asList("02c5c466ad26104b64166beb788a618c"));
        executeTest("test MultiSample Pilot2 with alleles passed in and emitting all sites", spec2);
    }

    @Test
    public void testSingleSamplePilot2() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,100,000", 1,
                Arrays.asList("7cc7ec249243a0dcb9e3a3033d4704e4"));
        executeTest("test SingleSample Pilot2", spec);
    }

    @Test
    public void testMultipleSNPAlleles() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " -nosl --no_cmdline_in_header -glm BOTH --dbsnp " + b37dbSNP129 + " -I " + testDir + "multiallelic.snps.bam -o %s -L " + testDir + "multiallelic.snps.intervals", 1,
                Arrays.asList("c0416efb6ff481889542e7807baa495c"));
        executeTest("test Multiple SNP alleles", spec);
    }

    @Test
    public void testBadRead() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " -nosl --no_cmdline_in_header -glm BOTH -I " + testDir + "badRead.test.bam -o %s -L 1:22753424-22753464", 1,
                Arrays.asList("e2cf97bca4a720ca64ca7f682da6c9f0"));
        executeTest("test bad read", spec);
    }

    @Test
    public void testReverseTrim() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " -nosl --no_cmdline_in_header -glm INDEL -I " + validationDataLocation + "CEUTrio.HiSeq.b37.chr20.10_11mb.bam -o %s -L 20:10289124 -L 20:10090289", 1,
                Arrays.asList("976694d7aedd0272a3e2ba46e00581f5"));
        executeTest("test reverse trim", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing compressed output
    //
    // --------------------------------------------------------------------------------------------------------------

    private final static String COMPRESSED_OUTPUT_MD5 = "c0f3d88386367fcd66d196b49a22bfaa";

    @Test
    public void testCompressedOutput() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,100,000", 1,
                Arrays.asList("gz"), Arrays.asList(COMPRESSED_OUTPUT_MD5));
        executeTest("test compressed output", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing parallelization
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testParallelization() {

        // Note that we need to turn off any randomization for this to work, so no downsampling and no annotations

        String md5 = "34cb7146c037925e8f324cffd986834d";

        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -dt NONE -G none -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,075,000", 1,
                Arrays.asList(md5));
        executeTest("test parallelization (single thread)", spec1);

        GenomeAnalysisEngine.resetRandomGenerator();

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -dt NONE -G none -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,075,000 -nt 2", 1,
                Arrays.asList(md5));
        executeTest("test parallelization (2 threads)", spec2);

        GenomeAnalysisEngine.resetRandomGenerator();

        WalkerTest.WalkerTestSpec spec3 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -dt NONE -G none -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,075,000 -nt 4", 1,
                Arrays.asList(md5));
        executeTest("test parallelization (4 threads)", spec3);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing parameters
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testMinBaseQualityScore() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 --min_base_quality_score 26", 1,
                Arrays.asList("ac4c8e0fdaf94a7949ceb157d3c836a9"));
        executeTest("test min_base_quality_score 26", spec);
    }

    @Test
    public void testSLOD() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b36KGReference + " --no_cmdline_in_header -glm BOTH --dbsnp " + b36dbSNP129 + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000", 1,
                Arrays.asList("751cfe684b0ed19b0d509bb7a315d65d"));
        executeTest("test SLOD", spec);
    }

    @Test
    public void testNDA() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " --annotateNDA -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000", 1,
                Arrays.asList("cbf1b8b261f726abe43045806e17e841"));
        executeTest("test NDA", spec);
    }

    @Test
    public void testCompTrack() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b36KGReference + " --no_cmdline_in_header -glm BOTH -comp:FOO " + b36dbSNP129 + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000", 1,
                Arrays.asList("70fd3b551f9def332393b7f3976440ac"));
        executeTest("test using comp track", spec);
    }

    @Test
    public void testOutputParameterSitesOnly() {
        testOutputParameters("-sites_only", "9700208609431548b2651dc582f1ea10");
    }

    @Test
    public void testOutputParameterAllConfident() {
        testOutputParameters("--output_mode EMIT_ALL_CONFIDENT_SITES", "b7bfc2a5298f86e3aef1bb70bebfb38e");
    }

    @Test
    public void testOutputParameterAllSites() {
        testOutputParameters("--output_mode EMIT_ALL_SITES", "90dcce8cd6de92258f332b261e64a5e9");
    }

    private void testOutputParameters(final String args, final String md5) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 " + args, 1,
                Arrays.asList(md5));
        executeTest(String.format("testParameter[%s]", args), spec);
    }

    @Test
    public void testConfidence() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 -stand_call_conf 10 ", 1,
                Arrays.asList("c5311739d1bba16cd59fb33e5f92ff49"));
        executeTest("test confidence 1", spec1);
    }

    @Test
    public void testConfidence2() {
        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 -stand_emit_conf 10 ", 1,
                Arrays.asList("189f5b7026595a0814e1eb1466834760"));
        executeTest("test confidence 2", spec2);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing heterozygosity
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testHeterozyosity1() {
        testHeterozosity( 0.01, "a818570e33316086ad5cb95992faca69" );
    }

    @Test
    public void testHeterozyosity2() {
        testHeterozosity( 1.0 / 1850, "3bfdba9bcc9c84d326750ad4a879897f" );
    }

    private void testHeterozosity(final double arg, final String md5) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,100,000 --heterozygosity " + arg, 1,
                Arrays.asList(md5));
        executeTest(String.format("test heterozyosity[%s]", arg), spec);
    }


    // --------------------------------------------------------------------------------------------------------------
    //
    // testing calls with SLX, 454, and SOLID data
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testMultiTechnologies() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand +
                        " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
                        " -o %s" +
                        " -L 1:10,000,000-10,100,000",
                1,
                Arrays.asList("1bd18653e8d9bf72c180142c2287eb5e"));

        executeTest(String.format("test multiple technologies"), spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing calls with BAQ
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testCallingWithBAQ() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand +
                        " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
                        " -o %s" +
                        " -L 1:10,000,000-10,100,000" +
                        " -baq CALCULATE_AS_NECESSARY",
                1,
                Arrays.asList("1defb30f6b870fbe5bcb27aed22ffb04"));

        executeTest(String.format("test calling with BAQ"), spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing indel caller
    //
    // --------------------------------------------------------------------------------------------------------------
    // Basic indel testing with SLX data
    @Test
    public void testSimpleIndels() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels +
                        " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam" +
                        " -o %s" +
                        " -L 1:10,000,000-10,500,000",
                1,
                Arrays.asList("2d52f5243da578c4070b03881852aae3"));

        executeTest(String.format("test indel caller in SLX"), spec);
    }

    // Basic indel testing with SLX data
    @Test
    public void testIndelsWithLowMinAlleleCnt() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels +
                        " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam" +
                        " -o %s" +
                        " -minIndelCnt 1" +
                        " -L 1:10,000,000-10,100,000",
                1,
                Arrays.asList("9271105e630ab39cf1c88b338da54594"));

        executeTest(String.format("test indel caller in SLX with low min allele count"), spec);
    }

    @Test
    public void testMultiTechnologyIndels() {
         WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                 baseCommandIndels +
                         " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
                         " -o %s" +
                         " -L 1:10,000,000-10,500,000",
                 1,
                 Arrays.asList("80f0221c47c30a03cb902b1c54c1eb12"));

         executeTest(String.format("test indel calling, multiple technologies"), spec);
     }

    @Test
    public void testWithIndelAllelesPassedIn1() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + testDir + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("ce40ff952025e806638156c3cd002927"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in", spec);
    }

    @Test
    public void testWithIndelAllelesPassedIn2() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles "
                        + testDir + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("4738fa9e68a010e05c64ae965ed94de3"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in and emitting all sites", spec);
    }

    @Test
    public void testMultiSampleIndels1() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10450700-10551000", 1,
                Arrays.asList("0b0114364680a23e125f469208b96c6f"));
        List<File> result = executeTest("test MultiSample Pilot1 CEU indels", spec1).getFirst();

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + result.get(0).getAbsolutePath() + " -I " + validationDataLocation +
                        "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10450700-10551000", 1,
                Arrays.asList("adcf53b8dcfde7f2c657745751549bfe"));
        executeTest("test MultiSample Pilot1 CEU indels using GENOTYPE_GIVEN_ALLELES", spec2);
    }

    @Test
    public void testGGAwithNoEvidenceInReads() {
        final String vcf = "small.indel.test.vcf";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndelsb37 + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -alleles " + testDir + vcf + " -I " + validationDataLocation +
                        "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam -o %s -L " + validationDataLocation + vcf, 1,
                Arrays.asList("3e3ac23846801c34acbf10a1a527264a"));
        executeTest("test GENOTYPE_GIVEN_ALLELES with no evidence in reads", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing SnpEff
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testSnpEffAnnotationRequestedWithoutRodBinding() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10,022,000-10,025,000 " +
                "-A SnpEff",
                1,
                UserException.class);
        executeTest("testSnpEffAnnotationRequestedWithoutRodBinding", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing SnpEff
    //
    // --------------------------------------------------------------------------------------------------------------

    final static String assessMinIndelFraction = baseCommandIndelsb37 + " -I " + validationDataLocation
            + "978604.bam -L 1:978,586-978,626 -o %s --sites_only -rf Sample -goodSM 7377 -goodSM 22-0022 -goodSM 134 -goodSM 344029-53 -goodSM 14030";
    
    @Test
    public void testMinIndelFraction0() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                assessMinIndelFraction + " -minIndelFrac 0.0", 1,
                Arrays.asList("eb9d8d1b19617190762f62bad139263a"));
        executeTest("test minIndelFraction 0.0", spec);
    }

    @Test
    public void testMinIndelFraction25() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                assessMinIndelFraction + " -minIndelFrac 0.25", 1,
                Arrays.asList("952f46506d3ef9c11cafcc072a00750a"));
        executeTest("test minIndelFraction 0.25", spec);
    }

    @Test
    public void testMinIndelFraction100() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                assessMinIndelFraction + " -minIndelFrac 1", 1,
                Arrays.asList("50a6774b7d8f71fe0e125c204d50ba84"));
        executeTest("test minIndelFraction 1.0", spec);
    }
}
