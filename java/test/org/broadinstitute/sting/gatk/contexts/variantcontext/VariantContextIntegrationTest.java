package org.broadinstitute.sting.gatk.contexts.variantcontext;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;
import java.util.List;
import java.io.File;

public class VariantContextIntegrationTest extends WalkerTest {
    private static String root = "-T TestVariantContext" +
            " -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
            " -B vcf,VCF,/humgen/gsa-hpprojects/GATK/data/Validation_Data/yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf" +
            " -R " + oneKGLocation + "reference/human_b36_both.fasta";

    static HashMap<String, String> expectations = new HashMap<String, String>();
    static {
        expectations.put("-L 1:1-10000 --printPerLocus", "39c035ae756eb176a7baffcd0f0fe0af");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly", "33ab797ac3de900e75fc0d9b01efe482");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsStartinAtCurrentPosition", "38619d704068072a4ccfd354652957a2");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly --onlyContextsStartinAtCurrentPosition", "bf64ab634126382813a6a6b29a5f47d8");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType SNP", "be087f53429974b4e505cd59f9363bfe");
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
    public void testLargeScaleConversion() {
        // this really just tests that we are seeing the same number of objects over all of chr1
        WalkerTestSpec spec = new WalkerTestSpec( root + " -L 1" + " -o %s",
                1, // just one output file
                Arrays.asList("3fcdd982df080e6abc0afaba6abdf386"));
         executeTest("testLargeScaleConversion", spec);
    }
}