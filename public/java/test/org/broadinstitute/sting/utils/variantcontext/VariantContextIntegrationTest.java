

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
            " -L 1:1-1,000,000 -B:dbsnp,vcf " + b36dbSNP129 +
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
        new VCITTest("--printPerLocus", "");
        new VCITTest("--printPerLocus --onlyContextsOfType SNP", "");
        new VCITTest("--printPerLocus --onlyContextsOfType INDEL", "");
        new VCITTest("--printPerLocus --onlyContextsOfType MIXED", "");
        new VCITTest("--printPerLocus --onlyContextsOfType NO_VARIATION", "");
        new VCITTest("--printPerLocus --takeFirstOnly", "");
        new VCITTest("--printPerLocus --onlyContextsOfType INDEL --onlyContextsStartinAtCurrentPosition", "");
        new VCITTest("--printPerLocus --onlyContextsStartinAtCurrentPosition", "");
        new VCITTest("--printPerLocus --takeFirstOnly --onlyContextsStartinAtCurrentPosition", "");

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
