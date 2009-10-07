package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SingleSampleGenotyperIntegrationTest extends WalkerTest {
    public static String baseTestString() {
        return "-T SingleSampleGenotyper -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout %s";
    }

    public static String testGeliLod5() {
        return baseTestString() + " --variant_output_format GELI -lod 5";
    }

    private static String OneMb1StateMD5 = "ca300adf2442d8874b4f13819547a7fb";
    private static String OneMb3StateMD5 = "45caee5f71a30c7175c7389b594369b1";
    private static String OneMbEmpiricalMD5 = "c8602b745b5c024f09938a9f87f923cf";

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
                "-T SingleSampleGenotyper" +
                        " -R /broad/1KG/reference/human_b36_both.fasta" +
                        " -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
                        " -varout %s" +
                        " -L 1:10,000,000-10,100,000" +
                        " -m empirical",
                1,
                Arrays.asList("09f7ea954a44e68b80726147aee10944"));

        executeTest(String.format("testMultiTechnologies"), spec);
    }    

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing the cache
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testCache() {
        for ( BaseMismatchModel model : BaseMismatchModel.values() ) {
            // calculated the expected value without the cache enabled
            WalkerTest.WalkerTestSpec withoutCacheSpec = new WalkerTest.WalkerTestSpec(
                    testGeliLod5() + " -L 1:10,000,000-10,100,000 --disableCache -m " + model.toString(), 1,
                    Arrays.asList(""));
            List<String> withoutCache = executeTest("empirical1MbTest", withoutCacheSpec ).getSecond();

            WalkerTest.WalkerTestSpec withCacheSpec = new WalkerTest.WalkerTestSpec(
                    testGeliLod5() + " -L 1:10,000,000-10,100,000 -m " + model.toString(), 1,
                    withoutCache);
            executeTest(String.format("testCache[%s]", model), withCacheSpec );
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing genotype mode
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void genotypeTest() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                testGeliLod5() + " -L 1:10,000,000-10,100,000 -m empirical --genotype", 1,
                Arrays.asList("b4c6d84ea14f9c34b04b799d3edd199a"));
        executeTest("genotypeTest", spec);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // basic base calling models
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void oneState100bpTest() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec( testGeliLod5() + " -L 1:10,000,000-10,000,100 -m one_state", 1, Arrays.asList("3cd402d889c015be4a318123468f4262"));
        executeTest("oneState100bpTest", spec);
    }

    @Test
    public void oneState1MbTest() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                testGeliLod5() + " -L 1:10,000,000-11,000,000 -m one_state",
                1, Arrays.asList(OneMb1StateMD5));
        executeTest("oneState1MbTest", spec);
    }

    @Test
    public void threeState1MbTest() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                testGeliLod5() + " -L 1:10,000,000-11,000,000 -m three_state", 1,
                Arrays.asList(OneMb3StateMD5));
        executeTest("threeState1MbTest", spec);
    }

    @Test
    public void empirical1MbTest() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                testGeliLod5() + " -L 1:10,000,000-11,000,000 -m empirical", 1,
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
        e.put( 10.0, "3c6c40551c42f693736dab5e29f7a76c" );
        e.put( 3.0, "5c82a81f5919b29058f9ca031c00b7ef" );

        for ( Map.Entry<Double, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    baseTestString() + " --variant_output_format GELI -L 1:10,000,000-11,000,000 -m EMPIRICAL -lod " + entry.getKey(), 1,
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
        e.put( 0.01, "c6fc697d31cfda563145138126561c77" );
        e.put( 0.0001, "1a1466b7599c60fad75d0513a0c8496e" );
        e.put( 1.0 / 1850, "5758dbcdea158fb053fcea22b06befe8" );
        
        for ( Map.Entry<Double, String> entry : e.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    testGeliLod5() + " -L 1:10,000,000-11,000,000 -m EMPIRICAL --heterozygosity " + entry.getKey(), 1,
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
