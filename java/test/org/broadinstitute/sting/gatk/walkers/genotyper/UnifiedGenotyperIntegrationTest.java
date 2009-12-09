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
    public static String baseTestString() {
        return "-T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s";
    }

    public static String testGeliLod5() {
        return baseTestString() + " --variant_output_format GELI -confidence 50";
    }

    private static String OneMb1StateMD5 = "7e3fc1d8427329eb2a3e05a81011749a";
    private static String OneMb3StateMD5 = "f5912d5d6585436a77495688f09cf1cc";
    private static String OneMbEmpiricalMD5 = "b9b2d9c7eb9a7af416fddd3b77a72efe";

//    private static String oneMbMD5(BaseMismatchModel m) {
//        switch (m) {
//            case ONE_STATE: return OneMb1StateMD5;
//            case THREE_STATE: return OneMb3StateMD5;
//            case EMPIRICAL: return OneMbEmpiricalMD5;
//            default: throw new RuntimeException("Unexpected BaseMismatchModel " + m);
//        }
//    }

    // Uncomment to not check outputs against expectations
    //protected boolean parameterize() {
    //    return true;
    //}

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing multi-sample calling
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testMultiSamplePilot1PointEM() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/low_coverage_CEU.chr1.10k-11k.bam -varout %s -L 1:10,023,400-10,024,000 -bm empirical -gm EM_POINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("12d07018d9604b86b854e0c9132b0608"));
        executeTest("testMultiSamplePilot1 - Point Estimate EM", spec);
    }

    @Test
    public void testMultiSamplePilot2PointEM() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/pilot2_daughters.chr20.10k-11k.bam -varout %s -L 20:10,000,000-10,010,000 -bm empirical -gm EM_POINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("65a808f077447f89e001e1cce8433acc"));
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
                "-T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/low_coverage_CEU.chr1.10k-11k.bam -varout %s -L 1:10,023,000-10,024,000 -bm empirical -gm POOLED -ps 60 -confidence 30", 1,
                Arrays.asList("4890ce37667379a0b60892cb0b9d8e09"));
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
                "-T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/low_coverage_CEU.chr1.10k-11k.bam -varout %s -L 1:10,022,000-10,025,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("fd5291c7115a558a7a6309ab417d2d57"));
        executeTest("testMultiSamplePilot1 - Joint Estimate", spec);
    }

    @Test
    public void testMultiSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/pilot2_daughters.chr20.10k-11k.bam -varout %s -L 20:10,000,000-10,050,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("58ddc9426a25d0c3604af8da4eaa8294"));
        executeTest("testMultiSamplePilot2 - Joint Estimate", spec);
    }

    @Test
    public void testSingleSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,100,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("52e7030fe0968cd1a9b95154d8056b3f"));
        executeTest("testSingleSamplePilot2 - Joint Estimate", spec);
    }

    @Test
    public void testGenotypeModeJoint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -genotype -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,001,000 -bm empirical -gm JOINT_ESTIMATE -confidence 70", 1,
                Arrays.asList("1ac9c9c8be0eb6adee7c67fac7c19c0d"));
        executeTest("testGenotypeMode - Joint Estimate", spec);
    }

    @Test
    public void testAllBasesModeJoint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -all_bases -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,001,000 -bm empirical -gm JOINT_ESTIMATE -confidence 70", 1,
                Arrays.asList("37177e194a602749d31e5ecbf80ad309"));
        executeTest("testAllBasesMode - Joint Estimate", spec);
    }

    //@Test
    //public void testGLF() {
    //    WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
    //            "-T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,050,000 -bm empirical -gm JOINT_ESTIMATE -confidence 10", 1,
    //            Arrays.asList("a95b871bc0bc984f66815b20db7467fe"));
    //    executeTest("testGLF", spec);
    //}

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing calls with SLX, 454, and SOLID data
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testMultiTechnologies() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper" +
                        " -R /broad/1KG/reference/human_b36_both.fasta" +
                        " -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
                        " -varout %s" +
                        " -L 1:10,000,000-10,100,000" +
                        " -bm empirical" +
		        " -vf GELI",
                1,
                Arrays.asList("46acd5a5d7301129ee3c0c9a4ab10259"));

        executeTest(String.format("testMultiTechnologies"), spec);
    }    

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing the cache
    //
    // --------------------------------------------------------------------------------------------------------------
    /*
    @Test
    public void testCache() {
        for ( BaseMismatchModel model : BaseMismatchModel.values() ) {
            // calculated the expected value without the cache enabled
            WalkerTest.WalkerTestSpec withoutCacheSpec = new WalkerTest.WalkerTestSpec(
                    testGeliLod5() + " -L 1:10,000,000-10,100,000 --disableCache -m " + model.toString(), 1,
                    Arrays.asList(""));
            List<String> withoutCache = executeTest("empirical1MbTest", withoutCacheSpec ).getSecond();

            WalkerTest.WalkerTestSpec withCacheSpec = new WalkerTest.WalkerTestSpec(
                    testGeliLod5() + " -L 1:10,000,000-10,100,000 -bm " + model.toString(), 1,
                    withoutCache);
            executeTest(String.format("testCache[%s]", model), withCacheSpec );
        }
    }
    */

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing genotype mode
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void genotypeTest() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                testGeliLod5() + " -L 1:10,000,000-10,100,000 -bm empirical --genotype", 1,
                Arrays.asList("db34a348abd13ee894780296e8a7e909"));
        executeTest("genotypeTest", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // basic base calling models
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void oneState100bpTest() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec( testGeliLod5() + " -L 1:10,000,000-10,000,100 -bm one_state", 1, Arrays.asList("3cd402d889c015be4a318123468f4262"));
        executeTest("oneState100bpTest", spec);
    }

    @Test
    public void oneState1MbTest() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                testGeliLod5() + " -L 1:10,000,000-11,000,000 -bm one_state",
                1, Arrays.asList(OneMb1StateMD5));
        executeTest("oneState1MbTest", spec);
    }

    @Test
    public void threeState1MbTest() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                testGeliLod5() + " -L 1:10,000,000-11,000,000 -bm three_state", 1,
                Arrays.asList(OneMb3StateMD5));
        executeTest("threeState1MbTest", spec);
    }

    @Test
    public void empirical1MbTest() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                testGeliLod5() + " -L 1:10,000,000-11,000,000 -bm empirical", 1,
                Arrays.asList(OneMbEmpiricalMD5));
        executeTest("empirical1MbTest", spec);
    }



    // --------------------------------------------------------------------------------------------------------------
    //
    // testing output formats
    //
    // --------------------------------------------------------------------------------------------------------------

    //@Argument(fullName = "variant_output_format", shortName = "vf", doc = "File format to be used", required = false)
    //public GenotypeWriterFactory.GENOTYPE_FORMAT VAR_FORMAT = GenotypeWriterFactory.GENOTYPE_FORMAT.GELI;

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing LOD thresholding
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testLOD() {
        HashMap<Double, String> e = new HashMap<Double, String>();
        e.put( 100.0, "94c5b48c0c956fcdacbffaa38a80d926" );
        e.put( 30.0, "df455abc1b2bc533aa1dc6eb088a835a" );

        for ( Map.Entry<Double, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    baseTestString() + " --variant_output_format GELI -L 1:10,000,000-11,000,000 -bm EMPIRICAL -confidence " + entry.getKey(), 1,
                    Arrays.asList(entry.getValue()));
            executeTest("testLOD", spec);
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing hetero setting
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testHeterozyosity() {
        HashMap<Double, String> e = new HashMap<Double, String>();
        e.put( 0.01, "b8837be7e8beb3ab2ed7150cdc022c65" );
        e.put( 0.0001, "ef0f2af7d13f166829d86b15fabc2b81" );
        e.put( 1.0 / 1850, "a435c8c966c11f4393a25a9d01c4fc3d" );

        for ( Map.Entry<Double, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    testGeliLod5() + " -L 1:10,000,000-11,000,000 -bm EMPIRICAL --heterozygosity " + entry.getKey(), 1,
                    Arrays.asList(entry.getValue()));
            executeTest(String.format("testHeterozyosity[%s]", entry.getKey()), spec);
        }
    }

    /**
     * test the output of a binary geli file
      */
    @Test
    public void empirical1MbTestBinaryGeli() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseTestString() + " -L 1:10,000,000-11,000,000 -bm empirical --variant_output_format GELI_BINARY -confidence 50", 1,
                Arrays.asList("18f175c7ccaeca57b8d412e9f4ebbe50"));
        executeTest("empirical1MbTestBinaryGeli", spec);
    }
}
