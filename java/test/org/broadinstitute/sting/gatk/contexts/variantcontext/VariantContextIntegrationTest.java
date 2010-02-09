package org.broadinstitute.sting.gatk.contexts.variantcontext;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;
import java.util.List;
import java.io.File;

public class VariantContextIntegrationTest extends WalkerTest {
    private static String cmdRoot = "-T TestVariantContext" +
            " -R " + oneKGLocation + "reference/human_b36_both.fasta";

    private static String root = cmdRoot +
            " -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
            " -B vcf,VCF,/humgen/gsa-hpprojects/GATK/data/Validation_Data/yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf";

    static HashMap<String, String> expectations = new HashMap<String, String>();
    static {
        expectations.put("-L 1:1-10000 --printPerLocus", "eb5802e7e615fa79b788a9447b7e824c");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly", "f4580866bcff10e0e742f422c45695d3");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsStartinAtCurrentPosition", "85c949c985d8aa340030362ea1fc13d2");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly --onlyContextsStartinAtCurrentPosition", "a45430cfed2f574fb69e79d74ed41017");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType SNP", "105547f432744e0776bd58e753cfa859");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType INDEL", "d758adbab9011e42c77d502fe4d62c27");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType INDEL --onlyContextsStartinAtCurrentPosition", "933ec8327192c5ed58a1952c73fb4f73");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType MIXED", "7d5d0283d92220ee78db7465d675b37f");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType NO_VARIATION", "39335acdb34c8a2af433dc50d619bcbc");
    }

    @Test
    public void testConversionSelection() {

        for ( Map.Entry<String, String> entry : expectations.entrySet() ) {
            String extraArgs = entry.getKey();
            String md5 = entry.getValue();

            WalkerTestSpec spec = new WalkerTestSpec( root + " " + extraArgs + " -o %s",
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testDbSNPAndVCFConversions", spec);
        }
    }

    @Test
    public void testToVCF() {
        // this really just tests that we are seeing the same number of objects over all of chr1

        WalkerTestSpec spec = new WalkerTestSpec( cmdRoot + " -B vcf,VCF," + validationDataLocation + "/yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.500.vcf -L 1:1-1000000 -o %s --outputVCF %s",
                2, // just one output file
                Arrays.asList("e3c35d0c4b5d4935c84a270f9df0951f", "461960b26ee1f8998ccc47da9bd3913c"));
         executeTest("testToVCF", spec);
    }

    @Test
    public void testLargeScaleConversion() {
        // this really just tests that we are seeing the same number of objects over all of chr1
        WalkerTestSpec spec = new WalkerTestSpec( root + " -L 1" + " -o %s",
                1, // just one output file
                Arrays.asList("3fcdd982df080e6abc0afaba6abdf386"));
         executeTest("testLargeScaleConversion", spec);
    }
}