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
                Arrays.asList("258e1954e6ae55c89abc6a716e19cbe0"));
        executeTest("test MultiSample Pilot1", spec);
    }

    //    @Test
    // todo - currently not working because when calling indels, using GENOTYPE_GIVEN_ALLELES yields a different result than in normal mode. To be fixed when extended events are removed.
    public void testMultiSamplePilot2AndRecallingWithAlleles() {
        String md5 = "b45636b29891f9df573ad2af6f507ee0";

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
                Arrays.asList("edeb1db288a24baff59575ceedd94243"));
        executeTest("test MultiSample Pilot2 with alleles passed in", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -B:alleles,vcf " + validationDataLocation + "allelesForUG.vcf -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,025,000", 1,
                Arrays.asList("581990130d90071b084024f4cd7caf91"));
        executeTest("test MultiSample Pilot2 with alleles passed in and emitting all sites", spec2);
    }

    @Test
    public void testSingleSamplePilot2() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,100,000", 1,
                Arrays.asList("d120db27d694a6da32367cc4fb5770fa"));
        executeTest("test SingleSample Pilot2", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing compressed output
    //
    // --------------------------------------------------------------------------------------------------------------

    private final static String COMPRESSED_OUTPUT_MD5 = "75e5c430ed39f79f24e375037a388dc4";

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

        String md5 = "a29615dd37222a11b8dadd341b53e43c";

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
        e.put( "--min_base_quality_score 26", "93e6269e38db9bc1732555e9969e3648" );
        e.put( "--min_mapping_quality_score 26", "64be99183c100caed4aa5f8bad64c7e9" );
        e.put( "--p_nonref_model GRID_SEARCH", "0592fe33f705ad8e2f13619fcf157805" );

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
        e.put( "-sites_only", "1483e637dc0279935a7f90d136d147bb" );
        e.put( "--output_mode EMIT_ALL_CONFIDENT_SITES", "adcd91bc7dae8020df8caf1a30060e98" );
        e.put( "--output_mode EMIT_ALL_SITES", "b708acc2fa40f336bcd2d0c70091e07e" );

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
                Arrays.asList("64be99183c100caed4aa5f8bad64c7e9"));
        executeTest("test confidence 1", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000 -stand_emit_conf 10 ", 1,
                Arrays.asList("e76ca54232d02f0d92730e1affeb804e"));
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
        e.put( 0.01, "18d37f7f107853b5e32c757b4e143205" );
        e.put( 1.0 / 1850, "2bcb90ce2f7542bf590f7612018fae8e" );

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
                Arrays.asList("825f05b31b5bb7e82231a15c7e4e2b0d"));

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
                Arrays.asList("0919ab7e513c377610e23a67d33608fa"));

        executeTest(String.format("test calling with BAQ"), spec);
    }

    @Test
    public void testCallingWithBAQOff() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand +
                        " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
                        " -o %s" +
                        " -L 1:10,000,000-10,100,000" +
                        " -baq OFF",
                1,
                Arrays.asList("825f05b31b5bb7e82231a15c7e4e2b0d"));

        executeTest(String.format("test calling with BAQ OFF"), spec);
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
                Arrays.asList("cb37348c41b8181be829912730f747e1"));

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
                Arrays.asList("ca5b6a5fb53ae401b146cc3044f454f2"));

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
                 Arrays.asList("ca4343a4ab6d3cce94ce61d7d1910f81"));

         executeTest(String.format("test indel calling, multiple technologies"), spec);
     }

    @Test
    public void testWithIndelAllelesPassedIn() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --genotyping_mode GENOTYPE_GIVEN_ALLELES -B:alleles,vcf " + validationDataLocation + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("3f555b53e9dd14cf7cdf96c24e322364"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                baseCommandIndels + " --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -B:alleles,vcf "
                        + validationDataLocation + "indelAllelesForUG.vcf -I " + validationDataLocation +
                        "pilot2_daughters.chr20.10k-11k.bam -o %s -L 20:10,000,000-10,100,000", 1,
                Arrays.asList("1b9764b783acf7822edc58e6822eef5b"));
        executeTest("test MultiSample Pilot2 indels with alleles passed in and emitting all sites", spec2);
    }


}
