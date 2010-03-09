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
    @Test
    public void testPooled1() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam -varout %s -L 1:10,023,000-10,024,000 -bm empirical -gm POOLED -ps 60 -confidence 30", 1,
                Arrays.asList("c30af5d192661abd77b05a316f1d8923"));
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
                Arrays.asList("882b2fae1cd1ba65cac3cadacec0ce2b"));
        executeTest("testMultiSamplePilot1 - Joint Estimate", spec);
    }

    @Test
    public void testMultiSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "pilot2_daughters.chr20.10k-11k.bam -varout %s -L 20:10,000,000-10,050,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("aa0cff414e6623c36465726a987a645d"));
        executeTest("testMultiSamplePilot2 - Joint Estimate", spec);
    }

    @Test
    public void testSingleSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,100,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("53df224164083cc7d8ad85f3d16ba38f"));
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
                Arrays.asList("3ade750c0d261594ea549db7b127a1e3"));
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
        e.put( "-genotype", "bee9fa71d70fdde094ab30785d4fa84e" );
        e.put( "-all_bases", "410cff9d97cd017becd1f6260c7abeeb" );
        e.put( "--min_base_quality_score 26", "85e1c35d3926afc68761aefea3f41332" );
        e.put( "--min_mapping_quality_score 26", "1c49a7d5e6ad295c0450b8a35053050f" );
        e.put( "--max_mismatches_in_40bp_window 5", "7e7db5a0d859704e12a4b89d35065682" );

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
                Arrays.asList("c67dd3e97cb188b117074d2c4692fcfa"));
        executeTest("testConfidence", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing beagle output
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testOtherOutput() {
        String[] md5s = {"ce0024816a092af9f998a7561ffb4fb2", "8cba0b8752f18fc620b4697840bc7291"};
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper" +
                        " -R " + oneKGLocation + "reference/human_b36_both.fasta" +
                        " -I " + validationDataLocation + "low_coverage_CEU.chr1.10k-11k.bam" +
                        " -varout %s" +
                        " -beagle %s" +
                        " -L 1:10,023,400-10,024,000" +
                        " -bm empirical" +
                        " -gm JOINT_ESTIMATE" +
		        " -vf VCF",
                2,
                Arrays.asList(md5s));

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
}
