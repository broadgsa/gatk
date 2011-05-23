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

    private final static String baseCommand = "-T UnifiedGenotyper -R " + b36KGReference + " -NO_HEADER -glm BOTH";
    private final static String baseCommandIndels = "-T UnifiedGenotyper -R " + b36KGReference + " -NO_HEADER -glm INDEL";

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing normal calling
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testMultiSamplePilot1() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -o %s -L 1:10,022,000-10,025,000", 1,
                Arrays.asList("14a3814607d8ff18ef40f49c454739f2"));
        executeTest("test MultiSample Pilot1", spec);
    }

    @Test
    public void testMultiSamplePilot2AndRecallingWithAlleles() {
        String md5 = "6aa91d95c360693e5c5f459c132bda32";

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
                Arrays.asList("8af463870c0f66cd3fccc5734ef86cb0"));
        executeTest("test MultiSample Pilot2 with alleles passed in", spec2);
    }

    @Test
    public void testSingleSamplePilot2() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,100,000", 1,
                Arrays.asList("6276a53a4dbf87c59dc84b356cca332f"));
        executeTest("test SingleSample Pilot2", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing compressed output
    //
    // --------------------------------------------------------------------------------------------------------------

    private final static String COMPRESSED_OUTPUT_MD5 = "af45d5a45ad2f030d777926717fc9049";

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
        e.put( "--min_base_quality_score 26", "93ad1ed265be26e0c1fe05975d01edd1" );
        e.put( "--min_mapping_quality_score 26", "01882d04939a5a6085fc5cadec019d3e" );
        e.put( "--p_nonref_model GRID_SEARCH", "7b0b7f536c66afe006174ff9a094dd95" );

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
        e.put( "-sites_only", "919d4b034ef0e07c8d50ea05b8fa0725" );
        e.put( "--output_mode EMIT_ALL_CONFIDENT_SITES", "a9bf545b6078eb772b9d01db4eaa4deb" );
        e.put( "--output_mode EMIT_ALL_SITES", "ab32c0f15401809d861cd8a0edb6f8ff" );

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
                Arrays.asList("01882d04939a5a6085fc5cadec019d3e"));
        executeTest("test confidence 1", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 -stand_emit_conf 10 ", 1,
                Arrays.asList("00807e8a5bc9172c34a5b4bd2e8afba7"));
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
        e.put( 0.01, "7db98e357bcebec5d362aacb0b176ca8" );
        e.put( 1.0 / 1850, "f3a25877f5ec5e0fa436a1d29b6a9b04" );

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
                Arrays.asList("b58d99ce5a71ed1dbd351cb509a83ac3"));

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
                Arrays.asList("032e33e33ede77d87c9442fb9895a2f8"));

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
                Arrays.asList("3a6ba2d9b9a5c606389d3353bb27bbe8"));

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
                Arrays.asList("f536ca23cb86b14fc20e453859dafc01"));

        executeTest(String.format("test indel caller in SLX witn low min allele count"), spec);
    }

    @Test
    public void testMultiTechnologyIndels() {
         WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                 baseCommandIndels +
                         " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
                         " -o %s" +
                         " -L 1:10,000,000-10,500,000",
                 1,
                 Arrays.asList("592a214b8ca1b62733f9627adb631f16"));

         executeTest(String.format("test indel calling, multiple technologies"), spec);
     }

    // todo - feature not yet fully working with indels
    //@Test
    public void testWithIndelAllelesPassedIn() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -B:alleles,vcf " + validationDataLocation + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("e95c545b8ae06f0721f260125cfbe1f0"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -B:alleles,vcf "
                        + validationDataLocation + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("6c96d76b9bc3aade0c768d7c657ae210"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in", spec2);
    }


}
