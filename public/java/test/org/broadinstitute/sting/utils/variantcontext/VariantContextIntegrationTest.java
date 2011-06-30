

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
        expectations.put("-L 1:1-10000 --printPerLocus", "e4ee2eaa3114888e918a1c82df7a027a");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly", "5b5635e4877d82e8a27d70dac24bda2f");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsStartinAtCurrentPosition", "ceced3f270b4fe407ee83bc9028becde");
        expectations.put("-L 1:1-10000 --printPerLocus --takeFirstOnly --onlyContextsStartinAtCurrentPosition", "9a9b9e283553c28bf58de1cafa38fe92");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType SNP", "2097e32988d603d3b353b50218c86d3b");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType INDEL", "033bd952fca048fe1a4f6422b57ab2ed");
        expectations.put("-L 1:1-10000 --printPerLocus --onlyContextsOfType INDEL --onlyContextsStartinAtCurrentPosition", "5e40980c02797f90821317874426a87a");
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
                Arrays.asList("529f936aa6c303658b23caf4e527782f"));
         executeTest("testLargeScaleConversion", spec);
    }
}
