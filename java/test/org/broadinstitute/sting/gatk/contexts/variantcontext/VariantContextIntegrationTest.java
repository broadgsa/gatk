

package org.broadinstitute.sting.gatk.contexts.variantcontext;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

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
        expectations.put("-L 1:1-10000 --printPerLocus", "fcc9b7ea66c4407f60317112c7d17aa0");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly", "897de66b05e10c80de01492c03842083");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsStartinAtCurrentPosition", "566ea10d38e8b7685ec95e4774a6fa05");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly --onlyContextsStartinAtCurrentPosition", "81fd480369c7e98a479a79b792d42305");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType SNP", "2097e32988d603d3b353b50218c86d3b");
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

        WalkerTestSpec spec = new WalkerTestSpec( cmdRoot + " -NO_HEADER -B:vcf,VCF " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.500.vcf -L 1:1-1000000 -o %s --outputVCF %s",
                2, // just one output file
                Arrays.asList("e3c35d0c4b5d4935c84a270f9df0951f", "e6673737acbb6bfabfcd92c4b2268241"));
         executeTest("testToVCF", spec);
    }

    @Test
    public void testLargeScaleConversion() {
        // this really just tests that we are seeing the same number of objects over all of chr1
        WalkerTestSpec spec = new WalkerTestSpec( root + " -L 1" + " -o %s",
                1, // just one output file
                Arrays.asList("d2a3f2fe329a0a64145cfd19fde45b99"));
         executeTest("testLargeScaleConversion", spec);
    }
}
