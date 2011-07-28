

package org.broadinstitute.sting.utils.variantcontext;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;

public class VariantContextIntegrationTest extends WalkerTest {
    private static String cmdRoot = "-T TestVariantContext" +
            " -R " + b36KGReference;

    private static String root = cmdRoot +
            " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
            " -B:vcf,VCF3 " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf";

    static HashMap<String, String> expectations = new HashMap<String, String>();
    static {
        expectations.put("-L 1:1-10000 --printPerLocus", "e9d96677a57bc3a10fb6d9ba942c19f0");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly", "8a1174d2b18b98e624abbe93e6af8fdd");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsStartinAtCurrentPosition", "3933f1fae5453c54c3f791a23de07599");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly --onlyContextsStartinAtCurrentPosition", "c9cf2f01bf045a58dcc7649fd6ea2396");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType SNP", "2097e32988d603d3b353b50218c86d3b");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType INDEL", "a103d856e8bc558c949c6e3f184e8913");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType INDEL --onlyContextsStartinAtCurrentPosition", "5f2265ac6c6d80d64dc6e69a05c1250b");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType MIXED", "06a3ae4c0afa23b429a9491ab7707f3c");
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

        WalkerTestSpec spec = new WalkerTestSpec( cmdRoot + " -NO_HEADER -B:vcf,VCF3 " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.500.vcf -L 1:1-1000000 -o %s --outputVCF %s",
                2, // just one output file
                Arrays.asList("e3c35d0c4b5d4935c84a270f9df0951f", "ff91731213fd0bbdc200ab6fd1c93e63"));
         executeTest("testToVCF", spec);
    }

    @Test
    public void testLargeScaleConversion() {
        // this really just tests that we are seeing the same number of objects over all of chr1
        WalkerTestSpec spec = new WalkerTestSpec( root + " -L 1" + " -o %s",
                1, // just one output file
                Arrays.asList("045a5b02c86aeb9301dc0b724da0c8f7"));
         executeTest("testLargeScaleConversion", spec);
    }
}
