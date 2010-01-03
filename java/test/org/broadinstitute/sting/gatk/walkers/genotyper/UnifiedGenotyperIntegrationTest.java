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
                Arrays.asList("8acbb0879a7c9f35e7cf7d77a7cab850"));
        executeTest("testMultiSamplePilot1 - Point Estimate EM", spec);
    }

    @Test
    public void testMultiSamplePilot2PointEM() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -varout %s -L 20:10,000,000-10,010,000 -bm empirical -gm EM_POINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("a5ec48e6695d9a5e7fcced7775baf6ff"));
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
                Arrays.asList("c20f370440ed55b0b6c38715985a850d"));
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
                Arrays.asList("0a2df3fc9fc72a3ace46aab899275636"));
        executeTest("testMultiSamplePilot1 - Joint Estimate", spec);
    }

    @Test
    public void testMultiSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -varout %s -L 20:10,000,000-10,050,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("52249b69d5d984a8bdc0438e1f106dcb"));
        executeTest("testMultiSamplePilot2 - Joint Estimate", spec);
    }

    @Test
    public void testSingleSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,100,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("172ca88b718f24fa0f2a3460b2dab33d"));
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
                Arrays.asList("d86d9148b1eee16da9a5bb57c929c0c9"));
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
        e.put( "-genotype", "2a8b8d4da5ceed494e2210b15c2c3eb0" );
        e.put( "-all_bases", "69b20ce7fc18f78011ee02b48a8d69aa" );
        e.put( "--min_base_quality_score 26", "67175e0262e626d46aaf3bc8680cf888" );
        e.put( "--min_mapping_quality_score 26", "8edb24206529a7a344008280ab69cbb2" );
        e.put( "--max_mismatches_in_40bp_window 5", "870812ed23b2c656439dc1c8e6376f89" );

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
                Arrays.asList("eb708a714a0ecbebed707f2d9684e165"));
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
        e.put( "GLF", "2dabd6313cbf65564933187511bf823e" );
        e.put( "GELI_BINARY", "e0af9777893c4b872b9ad867f0a2758c" );

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
        e.put( 0.01, "0556dc160a18eeaf1a7303d5c1e1ce46" );
        e.put( 1.0 / 1850, "0556dc160a18eeaf1a7303d5c1e1ce46" );

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
        e.put( "one_state", "78c4cfdd9ef21650709de5355178449a" );
        e.put( "three_state", "8ae5cebc06d520b1277b2f71265dd7e0" );

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
                Arrays.asList("1eabfec6a377992edb847d84608ab562"));

        executeTest(String.format("testMultiTechnologies"), spec);
    }
}
