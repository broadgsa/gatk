

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
            " -L 1:1-1,000,000 -B:dbsnp,vcf " + GATKDataLocation + "dbsnp_132.b36.excluding_sites_after_129.vcf" +
            " -B:vcf,VCF3 " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf";

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
        new VCITTest("--printPerLocus", "f36b81b8bcd210c0e3a1058d791b78ec");
        new VCITTest("--printPerLocus --onlyContextsOfType SNP", "a77492ba003a1fca8d8e0227fa642f34");
        new VCITTest("--printPerLocus --onlyContextsOfType INDEL", "9e0375a1b680d7df0971dbf256944d7a");
        new VCITTest("--printPerLocus --onlyContextsOfType MIXED", "93628cbba30033398e7e680b92cb3680");
        new VCITTest("--printPerLocus --onlyContextsOfType NO_VARIATION", "39335acdb34c8a2af433dc50d619bcbc");
        new VCITTest("--printPerLocus --takeFirstOnly", "c4a3d7545d26880635e0e5e4e69952e2");
        new VCITTest("--printPerLocus --onlyContextsOfType INDEL --onlyContextsStartinAtCurrentPosition", "22a7bb9e63d5f2950322c26397670e5c");
        new VCITTest("--printPerLocus --onlyContextsStartinAtCurrentPosition", "6387c1a400d1872ae4394d01e533c296");
        new VCITTest("--printPerLocus --takeFirstOnly --onlyContextsStartinAtCurrentPosition", "dde3a3db4d9c57f5042e0dfe03380987");

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

        WalkerTestSpec spec = new WalkerTestSpec( cmdRoot + " -NO_HEADER -B:vcf,VCF3 " + validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.500.vcf -L 1:1-1000000 -o %s --outputVCF %s",
                2, // just one output file
                Arrays.asList("e3c35d0c4b5d4935c84a270f9df0951f", "ff91731213fd0bbdc200ab6fd1c93e63"));
         executeTest("testToVCF", spec);
    }
}
