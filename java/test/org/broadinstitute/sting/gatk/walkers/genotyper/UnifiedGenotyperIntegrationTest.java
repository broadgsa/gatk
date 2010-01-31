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
                Arrays.asList("5ad3f97c886a3381821366caa9162c12"));
        executeTest("testMultiSamplePilot1 - Point Estimate EM", spec);
    }

    @Test
    public void testMultiSamplePilot2PointEM() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -varout %s -L 20:10,000,000-10,010,000 -bm empirical -gm EM_POINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("73c80566c8353958b2ac61932f0b3812"));
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
                Arrays.asList("5a86c53e8f0897a71ff74662c5696dc4"));
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
                Arrays.asList("5e0a92fbddeb9d6e35586d0488a1e5c7"));
        executeTest("testMultiSamplePilot1 - Joint Estimate", spec);
    }

    @Test
    public void testMultiSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -varout %s -L 20:10,000,000-10,050,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("68d00450e3c2129ea38c67171722b385"));
        executeTest("testMultiSamplePilot2 - Joint Estimate", spec);
    }

    @Test
    public void testSingleSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,100,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("8971ab1c9d2780e5e12e9bfc0b059cd1"));
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
                Arrays.asList("fb5d09eb8f1494d48032be7272699add"));
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
        e.put( "-genotype", "41af43f6eaab72de553d865a3089bf54" );
        e.put( "-all_bases", "7cc1609aef6d6cc3dd7822c52c403750" );
        e.put( "--min_base_quality_score 26", "9596e2102369ced181f2a87d686faf2e" );
        e.put( "--min_mapping_quality_score 26", "130efb2b8bd7b495bf65c6477bcf83c8" );
        e.put( "--max_mismatches_in_40bp_window 5", "18935308954cf390b628c9226eccbe94" );

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
                Arrays.asList("8c0e1fed37a2eac9eaaaa59e31350f43"));
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
                Arrays.asList("3c6d76d55d608482940cd725b87ef07d"));

        executeTest(String.format("testMultiTechnologies"), spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing beagle output
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testOtherOutput() {
        String[] md5s = {"78482125d51f9eb2ee850a6b25921e84", "cf39a4a90e9f01b6f1183e6b1fd7520e"};
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
