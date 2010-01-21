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
    // testing point estimate model
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testMultiSamplePilot1PointEM() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -varout %s -L 1:10,023,400-10,024,000 -bm empirical -gm EM_POINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("46232790dae2e99e79626de78836ba08"));
        executeTest("testMultiSamplePilot1 - Point Estimate EM", spec);
    }

    @Test
    public void testMultiSamplePilot2PointEM() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -varout %s -L 20:10,000,000-10,010,000 -bm empirical -gm EM_POINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("683cfa48cc5e2d140286645873184c20"));
        executeTest("testMultiSamplePilot2 - Point Estimate EM", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing pooled model
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testPooled1() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -varout %s -L 1:10,023,000-10,024,000 -bm empirical -gm POOLED -ps 60 -confidence 30", 1,
                Arrays.asList("f2b3799fe18010aa8cfd76b0e0782db7"));
        executeTest("testPooled1", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing joint estimation model
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testMultiSamplePilot1Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -varout %s -L 1:10,022,000-10,025,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("42fc8d585f2f30dc6c58d413738e1f7c"));
        executeTest("testMultiSamplePilot1 - Joint Estimate", spec);
    }

    @Test
    public void testMultiSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -varout %s -L 20:10,000,000-10,050,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("c854552bbf4a33b8dca488ea5f0a4a32"));
        executeTest("testMultiSamplePilot2 - Joint Estimate", spec);
    }

    @Test
    public void testSingleSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,100,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("9dcb1d4af7c02804150a0f6c38be4e1e"));
        executeTest("testSingleSamplePilot2 - Joint Estimate", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing joint estimation model
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testParallelization() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,400,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30 -nt 4", 1,
                Arrays.asList("35637b1ec52f2a7c0551b87795c3060c"));
        executeTest("test parallelization", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing parameters
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testParameter() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "-genotype", "e5d4287d27aa7734f0d57a8213f549ef" );
        e.put( "-all_bases", "3a291bb06764e3615230f467cc501096" );
        e.put( "--min_base_quality_score 26", "fb00499f249973c156ca64e30e7b2d91" );
        e.put( "--min_mapping_quality_score 26", "315c7951a783655445933ed9d83db2c2" );
        e.put( "--max_mismatches_in_40bp_window 5", "bca61f8afce9f64763a1225e9c614d38" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,010,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30 " + entry.getKey(), 1,
                    Arrays.asList(entry.getValue()));
            executeTest(String.format("testParameter[%s]", entry.getKey()), spec);
        }
    }

    @Test
    public void testConfidence() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,010,000 -bm empirical -gm JOINT_ESTIMATE -confidence 10 ", 1,
                Arrays.asList("c7427818f57cc3bf11e9ee98461c1a65"));
        executeTest("testConfidence", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing other output formats
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testOtherFormat() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "GLF", "d5d1c5ea8d42712d5509cd1c9c38359d" );
        e.put( "GELI_BINARY", "764a0fed1b3cf089230fd91f3be9c2df" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,100,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30 -vf " + entry.getKey(), 1,
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
                    "-T UnifiedGenotyper -vf GELI -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,100,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30 --heterozygosity " + entry.getKey(), 1,
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
                    "-T UnifiedGenotyper -vf GELI -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,100,000 -gm JOINT_ESTIMATE -confidence 30 -bm " + entry.getKey(), 1,
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
                        " -bm empirical" +
                        " -gm JOINT_ESTIMATE" +
		        " -vf GELI",
                1,
                Arrays.asList("f09ac61858c2633e5d1326fcf098b36d"));

        executeTest(String.format("testMultiTechnologies"), spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing beagle output
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testOtherOutput() {
        String[] md5s = {"a5dce541f00d3fe364d110f1cae53538", "0147670a1b62bb3d218d5ed3f9fc4656"};
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper" +
                        " -R " + oneKGLocation + "reference/human_b36_both.fasta" +
                        " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam" +
                        " -varout %s" +
                        " -beagle %s" +
                        " -L 1:10,023,400-10,024,000" +
                        " -bm empirical" +
                        " -gm JOINT_ESTIMATE" +
		        " -vf GELI",
                2,
                Arrays.asList(md5s));

        executeTest(String.format("testOtherOutput"), spec);
    }
}
