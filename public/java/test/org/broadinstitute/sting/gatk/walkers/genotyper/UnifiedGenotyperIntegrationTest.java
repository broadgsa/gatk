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
                Arrays.asList("cdec335abc9ad8e59335e39a73e0e95a"));
        executeTest("test MultiSample Pilot1", spec);
    }

    @Test
    public void testWithAllelesPassedIn1() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + privateTestDir + "allelesForUG.vcf -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,025,000", 1,
                Arrays.asList("efddb5e258f97fd4f6661cff9eaa57de"));
        executeTest("test MultiSample Pilot2 with alleles passed in", spec1);
    }

    @Test
    public void testWithAllelesPassedIn2() {
        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + privateTestDir + "allelesForUG.vcf -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,025,000", 1,
                Arrays.asList("24532eb381724cd74e99370da28d49ed"));
        executeTest("test MultiSample Pilot2 with alleles passed in and emitting all sites", spec2);
    }

    @Test
    public void testSingleSamplePilot2() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,100,000", 1,
                Arrays.asList("062a946160eec1d0fc135d58ca654ff4"));
        executeTest("test SingleSample Pilot2", spec);
    }

    @Test
    public void testMultipleSNPAlleles() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm BOTH --dbsnp " + b37dbSNP129 + " -I " + privateTestDir + "multiallelic.snps.bam -o %s -L " + privateTestDir + "multiallelic.snps.intervals", 1,
                Arrays.asList("94dc17d76d841f1d3a36160767ffa034"));
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
                Arrays.asList("9106d01ca0d0a8fedd068e72d509f380"));
        executeTest("test reverse trim", spec);
    }

    @Test
    public void testMismatchedPLs() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -glm INDEL -I " + privateTestDir + "mismatchedPLs.bam -o %s -L 1:24020341", 1,
                Arrays.asList("d847acf841ba8ba653f996ce4869f439"));
        executeTest("test mismatched PLs", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing compressed output
    //
    // --------------------------------------------------------------------------------------------------------------

    private final static String COMPRESSED_OUTPUT_MD5 = "6792419c482e767a3deb28913ed2b1ad";

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

        String md5 = "3d25ea660275a455ca443a786bff3d32";

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
                Arrays.asList("56157d930da6ccd224bce1ca93f11e41"));
        executeTest("test min_base_quality_score 26", spec);
    }

    @Test
    public void testSLOD() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b36KGReference + " --computeSLOD --no_cmdline_in_header -glm BOTH --dbsnp " + b36dbSNP129 + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000", 1,
                Arrays.asList("9f95cfe14d53a697c58247833bfd72a6"));
        executeTest("test SLOD", spec);
    }

    @Test
    public void testNDA() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " --annotateNDA -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000", 1,
                Arrays.asList("480437dd6e2760f4ab3194431519f331"));
        executeTest("test NDA", spec);
    }

    @Test
    public void testCompTrack() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b36KGReference + " --no_cmdline_in_header -glm BOTH -comp:FOO " + b36dbSNP129 + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000", 1,
                Arrays.asList("22c039412fd387dde6125b07c9a74a25"));
        executeTest("test using comp track", spec);
    }

    //@Test
    public void testNoCmdLineHeaderStdout() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandNoCmdLineHeaderStdout + " -glm INDEL -L 1:67,225,396-67,288,518", 0,
                Collections.<String>emptyList());
        executeTest("testNoCmdLineHeaderStdout", spec);
    }

    @Test
    public void testOutputParameterSitesOnly() {
        testOutputParameters("-sites_only", "40aeb4c9e31fe7046b72afc58e7599cb");
    }

    @Test
    public void testOutputParameterAllConfident() {
        testOutputParameters("--output_mode EMIT_ALL_CONFIDENT_SITES", "c706ca93b25ff83613cb4e95dcac567c");
    }

    @Test
    public void testOutputParameterAllSites() {
        testOutputParameters("--output_mode EMIT_ALL_SITES", "8a263fd0a94463ce1de9990f2b8ec841");
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
                Arrays.asList("df524e98903d96ab9353bee7c16a69de"));
        executeTest("test confidence 1", spec1);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing heterozygosity
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testHeterozyosity1() {
        testHeterozosity( 0.01, "8e61498ca03a8d805372a64c466b3b42" );
    }

    @Test
    public void testHeterozyosity2() {
        testHeterozosity( 1.0 / 1850, "668d06b5173cf3b97d052726988e1d7b" );
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
                Arrays.asList("908eb5e21fa39e7fb377cf4a9c4c7835"));

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
                Arrays.asList("c814558bb0ed2e19b12e1a2bf4465d52"));

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
                Arrays.asList("3593495aab5f6204c65de0b073a6ff65"));

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
                Arrays.asList("8b486a098029d5a106b0a37eff541c15"));

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
                 Arrays.asList("18efedc50cae2aacaba372265e38310b"));

         executeTest(String.format("test indel calling, multiple technologies"), spec);
     }

    @Test
    public void testWithIndelAllelesPassedIn1() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + privateTestDir + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("3ff8c7c80a518aa3eb8671a21479de5f"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in", spec);
    }

    @Test
    public void testWithIndelAllelesPassedIn2() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles "
                        + privateTestDir + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("578c0540f4f2052a634a829bcb9cc27d"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in and emitting all sites", spec);
    }

    @Test
    public void testMultiSampleIndels1() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10450700-10551000", 1,
                Arrays.asList("f7d0d0aee603df25c1f0525bb8df189e"));
        List<File> result = executeTest("test MultiSample Pilot1 CEU indels", spec1).getFirst();

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles " + result.get(0).getAbsolutePath() + " -I " + validationDataLocation +
                        "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10450700-10551000", 1,
                Arrays.asList("fc91d457a16b4ca994959c2b5f3f0352"));
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
                Arrays.asList("857b8e5df444463ac27f665c4f67fbe2"));
        executeTest("test minIndelFraction 0.0", spec);
    }

    @Test
    public void testMinIndelFraction25() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                assessMinIndelFraction + " -minIndelFrac 0.25", 1,
                Arrays.asList("81d4c7d9010fd6733b2997bc378e7471"));
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
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -I " + validationDataLocation + "testWithNs.bam -o %s -L 8:141799600-141814700", 1,
                Arrays.asList("bd7984a374f0ae5d277bd5fc5065f64f"));
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
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -I " + privateTestDir + "bamExample.ReducedRead.ADAnnotation.bam -o %s -L 1:67,225,396-67,288,518", 1,
                Arrays.asList("1b711c0966a2f76eb21801e73b76b758"));
        executeTest("test calling on a ReducedRead BAM", spec);
    }

    @Test
    public void testReducedBamSNPs() {
        testReducedCalling("SNP", "9ba4867cadb366746ee63e7a4afdb95e");
    }

    @Test
    public void testReducedBamINDELs() {
        testReducedCalling("INDEL", "132a4e0ccf9230b5bb4b56c649e2bdd5");
    }


    private void testReducedCalling(final String model, final String md5) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + b37KGReference + " --no_cmdline_in_header -I " + privateTestDir + "NA12878.HiSeq.b37.chr20.10_11mb.reduced.bam -o %s -L 20:10,000,000-11,000,000 -glm " + model, 1,
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
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 --contamination_fraction_to_filter 0.20", 1,
                Arrays.asList("27dd04159e06d9524fb8a4eef41f96ae"));
        executeTest("test contamination_percentage_to_filter 0.20", spec);
    }


}
