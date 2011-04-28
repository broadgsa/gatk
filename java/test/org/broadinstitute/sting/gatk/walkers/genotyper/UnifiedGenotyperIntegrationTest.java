package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

// ********************************************************************************** //
// Note that this class also serves as an integration test for the VariantAnnotator!  //
// ********************************************************************************** //

public class UnifiedGenotyperIntegrationTest extends WalkerTest {

    private final static String baseCommand = "-T UnifiedGenotyper -R " + b36KGReference + " -NO_HEADER";

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing normal calling
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testMultiSamplePilot1() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10,022,000-10,025,000", 1,
                Arrays.asList("a1add49e71306d72252ee3af3c40d2da"));
        executeTest("test MultiSample Pilot1", spec);
    }

    @Test
    public void testMultiSamplePilot2AndRecallingWithAlleles() {
        String md5 = "87a99063152ca935a1bec87ef19e0dad";

        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,050,000", 1,
                Arrays.asList(md5));
        List<File> result = executeTest("test MultiSample Pilot2", spec1).getFirst();

        GenomeAnalysisEngine.resetRandomGenerator();

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -B:alleles,vcf " + result.get(0).getAbsolutePath() + " -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,050,000", 1,
                Arrays.asList(md5));
        executeTest("test MultiSample Pilot2 with alleles passed in", spec2);       
    }

    @Test
    public void testWithAllelesPassedIn() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -B:alleles,vcf " + validationDataLocation + "allelesForUG.vcf -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,025,000", 1,
                Arrays.asList("a2584e8f92efa056cd8a09814476f883"));
        executeTest("test MultiSample Pilot2 with alleles passed in", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -B:alleles,vcf " + validationDataLocation + "allelesForUG.vcf -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,025,000", 1,
                Arrays.asList("0d15caf026c74bf859900e717b7fe996"));
        executeTest("test MultiSample Pilot2 with alleles passed in", spec2);
    }

    @Test
    public void testSingleSamplePilot2() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -glm SNP -L 1:10,000,000-10,100,000", 1,
                Arrays.asList("2241d97972b72233cb7b4c565739698e"));
        executeTest("test SingleSample Pilot2", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing compressed output
    //
    // --------------------------------------------------------------------------------------------------------------

    private final static String COMPRESSED_OUTPUT_MD5 = "64fd80527aa1b00fbbb4a77f1d5df48c";

    @Test
    public void testCompressedOutput() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,100,000", 1,
                Arrays.asList("gz"), Arrays.asList(COMPRESSED_OUTPUT_MD5));
        executeTest("test compressed output", spec);
    }

    // todo -- fixme
//    @Test
//    public void testCompressedOutputParallel() {
//        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
//                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,100,000 -nt 4", 1,
//                Arrays.asList("gz"), Arrays.asList(COMPRESSED_OUTPUT_MD5));
//        executeTest("testCompressedOutput-nt4", spec);
//    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing parallelization
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testParallelization() {

        // Note that we need to turn off any randomization for this to work, so no downsampling and no annotations

        String md5 = "f83a33a1ecc350cae0c002e4a43a7861";

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
    public void testCallingParameters() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "--min_base_quality_score 26", "f91280e3c12646207b3d06f9a563699b" );
        e.put( "--min_mapping_quality_score 26", "6f9e99469ec5f5ef333a5075e8149674" );
        e.put( "--p_nonref_model GRID_SEARCH", "430cfe737affe6b179725e1eb0b95655" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 " + entry.getKey(), 1,
                    Arrays.asList(entry.getValue()));
            executeTest(String.format("test calling parameter[%s]", entry.getKey()), spec);
        }
    }

    @Test
    public void testOutputParameter() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "-sites_only", "81e01f685e98d8f18fa288beb1246881" );
        e.put( "--output_mode EMIT_ALL_CONFIDENT_SITES", "3916157145aea448a38e4f12a0c381ae" );
        e.put( "--output_mode EMIT_ALL_SITES", "7652e2362e7fc038bc0f41a03df5d9f4" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 " + entry.getKey(), 1,
                    Arrays.asList(entry.getValue()));
            executeTest(String.format("testParameter[%s]", entry.getKey()), spec);
        }
    }

    @Test
    public void testConfidence() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 -stand_call_conf 10 ", 1,
                Arrays.asList("6f9e99469ec5f5ef333a5075e8149674"));
        executeTest("test confidence 1", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 -stand_emit_conf 10 ", 1,
                Arrays.asList("c2020de023bb3d731e17ed97aaa85c3d"));
        executeTest("test confidence 2", spec2);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing heterozygosity
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testHeterozyosity() {
        HashMap<Double, String> e = new HashMap<Double, String>();
        e.put( 0.01, "a56d3790bb2e18cec5a3099079431f19" );
        e.put( 1.0 / 1850, "119dff78b9f3b415ef3ed9d1bc46dd04" );

        for ( Map.Entry<Double, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,100,000 --heterozygosity " + entry.getKey(), 1,
                    Arrays.asList(entry.getValue()));
            executeTest(String.format("test heterozyosity[%s]", entry.getKey()), spec);
        }
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
                Arrays.asList("cd84c86ed02e4b052b96f76270ffd7fc"));

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
                Arrays.asList("9b71db6bebce80022c8fde5db09de65f"));

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
                baseCommand +
                        " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam" +
                        " -o %s" +
                        " -glm INDEL" +
                        " -L 1:10,000,000-10,500,000",
                1,
                Arrays.asList("84bc209f38d60f325f1a8b6292a82c82"));

        executeTest(String.format("test indel caller in SLX"), spec);
    }

    // Basic indel testing with SLX data
    @Test
    public void testIndelsWithLowMinAlleleCnt() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand +
                        " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam" +
                        " -o %s" +
                        " -glm INDEL -minIndelCnt 1" +
                        " -L 1:10,000,000-10,100,000",
                1,
                Arrays.asList("f536ca23cb86b14fc20e453859dafc01"));

        executeTest(String.format("test indel caller in SLX witn low min allele count"), spec);
    }

    @Test
    public void testMultiTechnologyIndels() {
         WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                 baseCommand +
                         " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
                         " -o %s" +
                         " -glm INDEL" +
                         " -L 1:10,000,000-10,500,000",
                 1,
                 Arrays.asList("e9b4fef2cdfa4c4657f0df53309131b6"));

         executeTest(String.format("test indel calling, multiple technologies"), spec);
     }

    // todo - feature not yet fully working with indels
    //@Test
    public void testWithIndelAllelesPassedIn() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -B:alleles,vcf " + validationDataLocation + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000 -glm INDEL", 1,
                Arrays.asList("e95c545b8ae06f0721f260125cfbe1f0"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -B:alleles,vcf "
                        + validationDataLocation + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000 -glm INDEL", 1,
                Arrays.asList("6c96d76b9bc3aade0c768d7c657ae210"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in", spec2);
    }


}
