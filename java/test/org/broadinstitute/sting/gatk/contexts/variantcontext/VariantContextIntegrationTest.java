

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
            " -B:vcf,VCF3 " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf";

    static HashMap<String, String> expectations = new HashMap<String, String>();
    static {
        expectations.put("-L 1:1-10000 --printPerLocus", "26ad6d6695a45c45e8451477fd9476a6");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly", "47772be91e4392d68aba901438aecdf2");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsStartinAtCurrentPosition", "e34688fa2450b0898ce55d0a5323db9a");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly --onlyContextsStartinAtCurrentPosition", "06505d1b90680862613ad374218b0d25");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType SNP", "2097e32988d603d3b353b50218c86d3b");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType INDEL", "8559d91b3f347c059d829fca1ada439e");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType INDEL --onlyContextsStartinAtCurrentPosition", "358330e2b373a38269abdf6e65180c0a");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType MIXED", "e5a00766f8c1ff9cf92310bafdec3126");
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
