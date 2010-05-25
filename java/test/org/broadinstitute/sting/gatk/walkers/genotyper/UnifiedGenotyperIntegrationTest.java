package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

// ********************************************************************************** //
// Note that this class also serves as an integration test for the VariantAnnotator!  //
// ********************************************************************************** //

public class UnifiedGenotyperIntegrationTest extends WalkerTest {

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing pooled model
    //
    // --------------------------------------------------------------------------------------------------------------
//    @Test
//    public void testPooled1() {
//        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
//                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -varout %s -L 1:10,023,000-10,024,000 -bm empirical -gm POOLED -ps 60 -confidence 30", 1,
//                Arrays.asList("c91f44a198cd7222520118726ea806ca"));
//        executeTest("testPooled1", spec);
//    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing joint estimation model
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testMultiSamplePilot1Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -varout %s -L 1:10,022,000-10,025,000", 1,
                Arrays.asList("3f919977584b2c15273cb7de948dd1ac"));
        executeTest("testMultiSamplePilot1 - Joint Estimate", spec);
    }

    @Test
    public void testMultiSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -varout %s -L 20:10,000,000-10,050,000", 1,
                Arrays.asList("d45934ae6537a84bb123db2288ad235e"));
        executeTest("testMultiSamplePilot2 - Joint Estimate", spec);
    }

    @Test
    public void testSingleSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,100,000", 1,
                Arrays.asList("96245f0842534527b66a9daad53134bd"));
        executeTest("testSingleSamplePilot2 - Joint Estimate", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing parallelization
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testParallelization() {
        String md5 = "638911a5d7b155076afc79a3d3f50548";

        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,075,000", 1,
                Arrays.asList(md5));
        executeTest("test parallelization (single thread)", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,075,000 -nt 2", 1,
                Arrays.asList(md5));
        executeTest("test parallelization (multithread)", spec2);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing parameters
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testParameter() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "-genotype", "d301500f0cd5ef1929a1adce7e857f9c" );
        e.put( "-all_bases", "70e94c6671a3e5122cb9ad04c0ad138e" );
        e.put( "--min_base_quality_score 26", "c508728f9d26cc3a8b73f3a5ee4275f2" );
        e.put( "--min_mapping_quality_score 26", "34634ea8f196eaf40dac6e069ebcf09a" );
        e.put( "--max_mismatches_in_40bp_window 5", "92903b97b43145c24f19a899ce127aa3" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,010,000 " + entry.getKey(), 1,
                    Arrays.asList(entry.getValue()));
            executeTest(String.format("testParameter[%s]", entry.getKey()), spec);
        }
    }

    @Test
    public void testConfidence() {
        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,010,000 -stand_call_conf 10 ", 1,
                Arrays.asList("cdcc3977ce78e15d7c8153eaf5f9a2df"));
        executeTest("testConfidence1", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,010,000 -stand_emit_conf 10 ", 1,
                Arrays.asList("2413aafdcfd3d8ac1169662f72524035"));
        executeTest("testConfidence2", spec2);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing beagle output
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testOtherOutput() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper" +
                        " -R " + oneKGLocation + "reference/human_b36_both.fasta" +
                        " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam" +
                        " -varout /dev/null" +
                        " -beagle %s" +
                        " -L 1:10,023,400-10,024,000",
                1,
                Arrays.asList("5077d80806e24865b5c12553843a5f85"));

        executeTest(String.format("testOtherOutput"), spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing other output formats
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testOtherFormat() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "GLF", "ddb1074b6f4a0fd1e15e4381476f1055" );
        e.put( "GELI_BINARY", "764a0fed1b3cf089230fd91f3be9c2df" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,100,000 -vf " + entry.getKey(), 1,
                    Arrays.asList(entry.getValue()));
            executeTest(String.format("testOtherFormat[%s]", entry.getKey()), spec);
        }
    }

    // --------------------------------------------- //
    // ALL REMAINING TESTS ARE OUTPUT IN GELI FORMAT //
    // --------------------------------------------- //

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing heterozygosity
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testHeterozyosity() {
        HashMap<Double, String> e = new HashMap<Double, String>();
        e.put( 0.01, "ee390f91867e8729b96220115e56ddb3" );
        e.put( 1.0 / 1850, "f96ad0ed71449bdb16b0c5561303a05a" );

        for ( Map.Entry<Double, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-T UnifiedGenotyper -vf GELI -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,100,000 --heterozygosity " + entry.getKey(), 1,
                    Arrays.asList(entry.getValue()));
            executeTest(String.format("testHeterozyosity[%s]", entry.getKey()), spec);
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing other base calling models
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testOtherBaseCallModel() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "one_state", "bcc983210b576d9fd228a67c5b9f372a" );
        e.put( "three_state", "2db3a5f3d46e13e2f44c34fbb7e7936f" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-T UnifiedGenotyper -vf GELI -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,100,000 -bm " + entry.getKey(), 1,
                    Arrays.asList(entry.getValue()));
            executeTest(String.format("testOtherBaseCallModel[%s]", entry.getKey()), spec);
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
                "-T UnifiedGenotyper" +
                        " -R " + oneKGLocation + "reference/human_b36_both.fasta" +
                        " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
                        " -varout %s" +
                        " -L 1:10,000,000-10,100,000" +
		                " -vf GELI",
                1,
                Arrays.asList("f67c690bf2e4eee2bb7c58b6646a2a98"));

        executeTest(String.format("testMultiTechnologies"), spec);
    }
}
