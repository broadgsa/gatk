

package org.broadinstitute.sting.utils.variantcontext;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;

public class VariantContextIntegrationTest extends WalkerTest {
    private static String cmdRoot = "-T TestVariantContext" +
            " -R " + b36KGReference;

    private static String root = cmdRoot +
            " -L 1:1-1,000,000 -V " + b36dbSNP129;

    private static final class VCITTest extends TestDataProvider {
        String args, md5;

        private VCITTest(final String args, final String md5) {
            super(VCITTest.class);
            this.args = args;
            this.md5 = md5;
        }
    }

    @DataProvider(name = "VCITTestData")
    public Object[][] createVCITTestData() {
        new VCITTest("--printPerLocus", "e9d0f1fe80659bb55b40aa6c3a2e921e");
        new VCITTest("--printPerLocus --onlyContextsOfType SNP", "0e620db3e45771df42c54a9c0ae4a29f");
        new VCITTest("--printPerLocus --onlyContextsOfType INDEL", "b725c204fefe3814644d50e7c20f9dfe");
        new VCITTest("--printPerLocus --onlyContextsOfType MIXED", "3ccc33f496a1718df55722d11cc14334");
        new VCITTest("--printPerLocus --onlyContextsOfType NO_VARIATION", "39335acdb34c8a2af433dc50d619bcbc");
        new VCITTest("--printPerLocus --takeFirstOnly", "3a45561da042b2b44b6a679744f16103");
        new VCITTest("--printPerLocus --onlyContextsOfType INDEL --onlyContextsStartinAtCurrentPosition", "4746f269ecc377103f83eb61cc162c39");
        new VCITTest("--printPerLocus --onlyContextsStartinAtCurrentPosition", "2749e3fae458650a85a2317e346dc44c");
        new VCITTest("--printPerLocus --takeFirstOnly --onlyContextsStartinAtCurrentPosition", "9bd48c2a40813023e29ffaa23d59d382");

        return VCITTest.getTests(VCITTest.class);
    }

    @Test(dataProvider = "VCITTestData")
    public void testConversionSelection(VCITTest test) {
	String extraArgs = test.args;
	String md5 = test.md5;

	WalkerTestSpec spec = new WalkerTestSpec( root + " " + extraArgs + " -o %s",
						  1, // just one output file
						  Arrays.asList(md5));
	executeTest("testSelectors", spec);
    }

    @Test
    public void testToVCF() {
        // this really just tests that we are seeing the same number of objects over all of chr1

        WalkerTestSpec spec = new WalkerTestSpec( cmdRoot + " -NO_HEADER -V:VCF3 " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.500.vcf -L 1:1-1000000 -o %s --outputVCF %s",
                2, // just one output file
                Arrays.asList("e3c35d0c4b5d4935c84a270f9df0951f", "ff91731213fd0bbdc200ab6fd1c93e63"));
         executeTest("testToVCF", spec);
    }
}
