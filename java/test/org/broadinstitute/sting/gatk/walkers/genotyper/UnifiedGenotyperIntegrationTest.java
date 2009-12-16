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

    private static String OneMb1StateMD5 = "c664bd887b89e9ebc0b4b569ad8eb128";
    private static String OneMb3StateMD5 = "9c68f6e900d081023ea97ec467a95bd8";
    private static String OneMbEmpiricalMD5 = "0b891eefeb2a2bfe707a8a0838b6d049";

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
                Arrays.asList("b19f85fddd08485c668e7192df79d944"));
        executeTest("testMultiSamplePilot1 - Point Estimate EM", spec);
    }

    @Test
    public void testMultiSamplePilot2PointEM() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/pilot2_daughters.chr20.10k-11k.bam -varout %s -L 20:10,000,000-10,010,000 -bm empirical -gm EM_POINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("a6c8e77d1741f3f6d958a0634fa59e14"));
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
                Arrays.asList("459b8f6a6265b67718495e9bffe96fa8"));
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
                Arrays.asList("c132c93cec5fdf02c3235180a7aa7dcc"));
        executeTest("testMultiSamplePilot1 - Joint Estimate", spec);
    }

    @Test
    public void testMultiSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/pilot2_daughters.chr20.10k-11k.bam -varout %s -L 20:10,000,000-10,050,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("09d48c56eae25fb6418e10feeb5b3fc5"));
        executeTest("testMultiSamplePilot2 - Joint Estimate", spec);
    }

    @Test
    public void testSingleSamplePilot2Joint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,100,000 -bm empirical -gm JOINT_ESTIMATE -confidence 30", 1,
                Arrays.asList("3544b1eb97834e2050239c101eaddb2d"));
        executeTest("testSingleSamplePilot2 - Joint Estimate", spec);
    }

    @Test
    public void testGenotypeModeJoint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -genotype -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,001,000 -bm empirical -gm JOINT_ESTIMATE -confidence 70", 1,
                Arrays.asList("7e8422f8008a4fe24ea1f5c2913b5d31"));
        executeTest("testGenotypeMode - Joint Estimate", spec);
    }

    @Test
    public void testAllBasesModeJoint() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T UnifiedGenotyper -all_bases -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s -L 1:10,000,000-10,001,000 -bm empirical -gm JOINT_ESTIMATE -confidence 70", 1,
                Arrays.asList("aedd59f9f2cae7fb07a53812b680852b"));
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
                Arrays.asList("eaca4b2323714dbd7c3ed379ce1843ba"));

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
                Arrays.asList("f9bdd9a8864467dbc4e5356bb8801a33"));
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
        e.put( 100.0, "6eec841b28fae433015b3d85608e03f7" );
        e.put( 30.0, "1b3365f41bbf6867516699afe9efc5f8" );

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
        e.put( 0.01, "601c48fc350083d14534ba5c3093edb9" );
        e.put( 0.0001, "bd03f7307314e45951d4d3e85fe63d16" );
        e.put( 1.0 / 1850, "662d479f1cd54480da1d0e66c81259b0" );

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
                Arrays.asList("b1027cf309c9ab7572528ce986e2c2d4"));
        executeTest("empirical1MbTestBinaryGeli", spec);
    }
}
