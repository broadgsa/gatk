package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class SingleSampleGenotyperIntegrationTest extends WalkerTest {
    public static String baseTestString() {
        return "-T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s";
    }

    public static String testGeliLod5() {
        return baseTestString() + " --variant_output_format GELI -lod 5";
    }

    private static String OneMb1StateMD5 = "cba4436066b5d88a461cc3eba74ad944";
    private static String OneMb3StateMD5 = "799dde914dd4afbffca9502b7f284780";
    private static String OneMbEmpiricalMD5 = "2b99b3c667c8ac94e8268b17e6979073";

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
                        " -bm empirical",
                1,
                Arrays.asList("0a923f676111b9bc27ccb1a9b97eafcc"));

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
                Arrays.asList("436c0e3365f61bf1d06eb630c025e51b"));
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
        e.put( 10.0, "3ec7815482112e0b9b4487ec69a52b67" );
        e.put( 3.0, "20da42fdce2e1c97d9c8d4935d31125d" );

        for ( Map.Entry<Double, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    baseTestString() + " --variant_output_format GELI -L 1:10,000,000-11,000,000 -bm EMPIRICAL -lod " + entry.getKey(), 1,
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
        e.put( 0.01, "e7494e6a62bd91ca02537c327a104395" );
        e.put( 0.0001, "1d51c711353d0db5b54dd0d2a7899c49" );
        e.put( 1.0 / 1850, "873f558a2c07ec40635e35d275b12d69" );

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
                baseTestString() + " -L 1:10,000,000-11,000,000 -m empirical --variant_output_format GELI_BINARY -lod 5", 1,
                Arrays.asList("17e9fc04b2a05cafc53562997c28e127"));
        executeTest("empirical1MbTest", spec);
    }

}
