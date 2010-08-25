
package org.broadinstitute.sting.gatk.contexts.variantcontext;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;

public class VariantContextIntegrationTest extends WalkerTest {
    private static String cmdRoot = "-T TestVariantContext" +
            " -R " + b36KGReference;

    private static String root = cmdRoot +
            " -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
            " -B:vcf,VCF " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf";

    static HashMap<String, String> expectations = new HashMap<String, String>();
    static {
        expectations.put("-L 1:1-10000 --printPerLocus", "63fd69e4ab430b79fb213dd27b58ae1c");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly", "276ed96efaaffc2fc1c3b3deb4e04d1d");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsStartinAtCurrentPosition", "a37f7bc34c1824688d3e475945c19d5a");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly --onlyContextsStartinAtCurrentPosition", "1715a6e0daf873f2e2cd10cb56085174");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType SNP", "bf33ab1ed65da7f56c02ca7956d9c31e");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType INDEL", "629ffd6b3b9ea1bce29cb715576f5c8a");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType INDEL --onlyContextsStartinAtCurrentPosition", "d4b812b2fec231f8f5b61d6f26cf86a5");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType MIXED", "546e8e546f2cdfba31f91ed083137c42");
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

        WalkerTestSpec spec = new WalkerTestSpec( cmdRoot + " -B:vcf,VCF " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.500.vcf -L 1:1-1000000 -o %s --outputVCF %s",
                2, // just one output file
                Arrays.asList("e3c35d0c4b5d4935c84a270f9df0951f", "f4db5f7346792b1155693722bc190f63"));
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
