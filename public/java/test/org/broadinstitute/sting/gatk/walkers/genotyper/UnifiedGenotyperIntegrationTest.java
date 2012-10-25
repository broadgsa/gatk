package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

// ********************************************************************************** //
// Note that this class also serves as an integration test for the VariantAnnotator!  //
// ********************************************************************************** //

public class UnifiedGenotyperIntegrationTest extends WalkerTest {

    private final static String baseCommand = "-T UnifiedGenotyper -R " + b36KGReference + " --no_cmdline_in_header -glm BOTH -minIndelFrac 0.0 --dbsnp " + b36dbSNP129;
    private final static String baseCommandIndels = "-T UnifiedGenotyper -R " + b36KGReference + " --no_cmdline_in_header -glm INDEL -mbq 20 -minIndelFrac 0.0 --dbsnp " + b36dbSNP129;
    private final static String baseCommandIndelsb37 = "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm INDEL -mbq 20 --dbsnp " + b37dbSNP132;
    private final static String baseCommandNoCmdLineHeaderStdout = "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -I " + privateTestDir + "bamExample.ReducedRead.ADAnnotation.bam";

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing normal calling
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testMultiSamplePilot1() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10,022,000-10,025,000", 1,
                Arrays.asList("847605f4efafef89529fe0e496315edd"));
        executeTest("test MultiSample Pilot1", spec);
    }

    @Test
    public void testWithAllelesPassedIn1() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + privateTestDir + "allelesForUG.vcf -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,025,000", 1,
                Arrays.asList("bc15123620e1134f799005d71d1180fe"));
        executeTest("test MultiSample Pilot2 with alleles passed in", spec1);
    }

    @Test
    public void testWithAllelesPassedIn2() {
        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + privateTestDir + "allelesForUG.vcf -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,025,000", 1,
                Arrays.asList("1ba7afccc8552f20d72d0b62237558e3"));
        executeTest("test MultiSample Pilot2 with alleles passed in and emitting all sites", spec2);
    }

    @Test
    public void testSingleSamplePilot2() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,100,000", 1,
                Arrays.asList("afb8768f31ab57eb43f75c1115eadc99"));
        executeTest("test SingleSample Pilot2", spec);
    }

    @Test
    public void testMultipleSNPAlleles() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm BOTH --dbsnp " + b37dbSNP129 + " -I " + privateTestDir + "multiallelic.snps.bam -o %s -L " + privateTestDir + "multiallelic.snps.intervals", 1,
                Arrays.asList("543f68e42034bf44cfb24da8c9204320"));
        executeTest("test Multiple SNP alleles", spec);
    }

    @Test
    public void testBadRead() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm BOTH -I " + privateTestDir + "badRead.test.bam -o %s -L 1:22753424-22753464", 1,
                Arrays.asList("d915535c1458733f09f82670092fcab6"));
        executeTest("test bad read", spec);
    }

    @Test
    public void testReverseTrim() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm INDEL -I " + validationDataLocation + "CEUTrio.HiSeq.b37.chr20.10_11mb.bam -o %s -L 20:10289124 -L 20:10090289", 1,
                Arrays.asList("5ce03dd9ca2d9324c1d4a9d64389beb5"));
        executeTest("test reverse trim", spec);
    }

    @Test
    public void testMismatchedPLs() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm INDEL -I " + privateTestDir + "mismatchedPLs.bam -o %s -L 1:24020341", 1,
                Arrays.asList("3c006b06b17bbe8e787d64eff6a63a19"));
        executeTest("test mismatched PLs", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing compressed output
    //
    // --------------------------------------------------------------------------------------------------------------

    private final static String COMPRESSED_OUTPUT_MD5 = "fd236bd635d514e4214d364f45ec4d10";

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

        String md5 = "d408b4661b820ed86272415b8ea08780";

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
                Arrays.asList("839ecd30d354a36b5dfa2b5e99859765"));
        executeTest("test min_base_quality_score 26", spec);
    }

    @Test
    public void testSLOD() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b36KGReference + " --computeSLOD --no_cmdline_in_header -glm BOTH --dbsnp " + b36dbSNP129 + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000", 1,
                Arrays.asList("c7429e670ba477bf9a6bbee2fb41c5a9"));
        executeTest("test SLOD", spec);
    }

    @Test
    public void testNDA() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " --annotateNDA -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000", 1,
                Arrays.asList("abd8e33e649cc11b55e200d3940cc7e2"));
        executeTest("test NDA", spec);
    }

    @Test
    public void testCompTrack() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b36KGReference + " --no_cmdline_in_header -glm BOTH -comp:FOO " + b36dbSNP129 + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000", 1,
                Arrays.asList("a6de72c5197f152be98a6f9029091a9e"));
        executeTest("test using comp track", spec);
    }

    @Test
    public void testNoCmdLineHeaderStdout() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandNoCmdLineHeaderStdout + " -glm INDEL -L 1:67,225,396-67,288,518", 0,
                Collections.<String>emptyList());
        executeTest("testNoCmdLineHeaderStdout", spec);
    }

    @Test
    public void testOutputParameterSitesOnly() {
        testOutputParameters("-sites_only", "97ba874eafc9884a4de027a84c036311");
    }

    @Test
    public void testOutputParameterAllConfident() {
        testOutputParameters("--output_mode EMIT_ALL_CONFIDENT_SITES", "f9ea04d96eeef29e71d37e60518c2579");
    }

    @Test
    public void testOutputParameterAllSites() {
        testOutputParameters("--output_mode EMIT_ALL_SITES", "41c046d38ea328421df924e37e017645");
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
                Arrays.asList("9addd225a985178339a0c49dc5fdc220"));
        executeTest("test confidence 1", spec1);
    }

    @Test
    public void testConfidence2() {
        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 -stand_emit_conf 10 ", 1,
                Arrays.asList("9addd225a985178339a0c49dc5fdc220"));
        executeTest("test confidence 2", spec2);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing heterozygosity
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testHeterozyosity1() {
        testHeterozosity( 0.01, "986923de51c71635d47e3d06fe3794a1" );
    }

    @Test
    public void testHeterozyosity2() {
        testHeterozosity( 1.0 / 1850, "fb12b1553f813004a394a391a8540873" );
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
                Arrays.asList("04a87b87ee4323eba853c78f25551d1a"));

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
                Arrays.asList("98058fc913b61c22d44875da1f5ea89c"));

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
                Arrays.asList("650c53774afacfc07a595675e8cdde17"));

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
                Arrays.asList("6a0c2a3a7bcc56ad01428c71408055aa"));

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
                 Arrays.asList("5f2721c3323de5390d2d47446139f32b"));

         executeTest(String.format("test indel calling, multiple technologies"), spec);
     }

    @Test
    public void testWithIndelAllelesPassedIn1() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + privateTestDir + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("7e3f67bf371112be5dbadb4fe6faa52a"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in", spec);
    }

    @Test
    public void testWithIndelAllelesPassedIn2() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles "
                        + privateTestDir + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("bc31c4977cb7e563ddf9c8dea27f3f4f"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in and emitting all sites", spec);
    }

    @Test
    public void testMultiSampleIndels1() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10450700-10551000", 1,
                Arrays.asList("a4761d7f25e7a62f34494801c98a0da7"));
        List<File> result = executeTest("test MultiSample Pilot1 CEU indels", spec1).getFirst();

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + result.get(0).getAbsolutePath() + " -I " + validationDataLocation +
                        "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10450700-10551000", 1,
                Arrays.asList("c526c234947482d1cd2ffc5102083a08"));
        executeTest("test MultiSample Pilot1 CEU indels using GENOTYPE_GIVEN_ALLELES", spec2);
    }

    @Test
    public void testGGAwithNoEvidenceInReads() {
        final String vcf = "small.indel.test.vcf";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndelsb37 + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -alleles " + privateTestDir + vcf + " -I " + validationDataLocation +
                        "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam -o %s -L " + validationDataLocation + vcf, 1,
                Arrays.asList("d76eacc4021b78ccc0a9026162e814a7"));
        executeTest("test GENOTYPE_GIVEN_ALLELES with no evidence in reads", spec);
    }

    @Test
    public void testBaseIndelQualityScores() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndelsb37 +
                        " -I " + privateTestDir + "NA12878.100kb.BQSRv2.example.bam" +
                        " -o %s" +
                        " -L 20:10,000,000-10,100,000",
                1,
                Arrays.asList("1e0d2c15546c3b0959b00ffb75488b56"));

        executeTest(String.format("test UG with base indel quality scores"), spec);
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
    // testing MinIndelFraction
    //
    // --------------------------------------------------------------------------------------------------------------

    final static String assessMinIndelFraction = baseCommandIndelsb37 + " -I " + validationDataLocation
            + "978604.bam -L 1:978,586-978,626 -o %s --sites_only -rf Sample -goodSM 7377 -goodSM 22-0022 -goodSM 134 -goodSM 344029-53 -goodSM 14030";

    @Test
    public void testMinIndelFraction0() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                assessMinIndelFraction + " -minIndelFrac 0.0", 1,
                Arrays.asList("ba4fafec383fb988f20c8cf53dd3e9a0"));
        executeTest("test minIndelFraction 0.0", spec);
    }

    @Test
    public void testMinIndelFraction25() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                assessMinIndelFraction + " -minIndelFrac 0.25", 1,
                Arrays.asList("4c57a88de275105156aaafc6f9041365"));
        executeTest("test minIndelFraction 0.25", spec);
    }

    @Test
    public void testMinIndelFraction100() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                assessMinIndelFraction + " -minIndelFrac 1", 1,
                Arrays.asList("3f07efb768e08650a7ce333edd4f9a52"));
        executeTest("test minIndelFraction 1.0", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing Ns in CIGAR
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testNsInCigar() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " -nosl --no_cmdline_in_header -I " + validationDataLocation + "testWithNs.bam -o %s -L 8:141799600-141814700", 1,
                Arrays.asList("e8ebfaac0804b782f22ab8ea35152735"));
        executeTest("test calling on reads with Ns in CIGAR", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing reduced reads
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testReducedBam() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " -nosl --no_cmdline_in_header -I " + privateTestDir + "bamExample.ReducedRead.ADAnnotation.bam -o %s -L 1:67,225,396-67,288,518", 1,
                Arrays.asList("da9c05f87bd6415e97f90c49cf68ed19"));
        executeTest("test calling on a ReducedRead BAM", spec);
    }

    @Test
    public void testReducedBamSNPs() {
        testReducedCalling("SNP", "1d4a826b144723ff0766c36aa0239287");
    }

    @Test
    public void testReducedBamINDELs() {
        testReducedCalling("INDEL", "68ef51d5c98480e0c0192e0eecb95bca");
    }


    private void testReducedCalling(final String model, final String md5) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " -nosl --no_cmdline_in_header -I " + privateTestDir + "NA12878.HiSeq.b37.chr20.10_11mb.reduced.bam -o %s -L 20:10,000,000-11,000,000 -glm " + model, 1,
                Arrays.asList(md5));
        executeTest("test calling on a ReducedRead BAM with " + model, spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing contamination down-sampling
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testContaminationDownsampling() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 --contamination_percentage_to_filter 0.20", 1,
                Arrays.asList("27dd04159e06d9524fb8a4eef41f96ae"));
        executeTest("test contamination_percentage_to_filter 0.20", spec);
    }


}
