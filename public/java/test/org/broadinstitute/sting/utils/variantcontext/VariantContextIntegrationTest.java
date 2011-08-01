

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
            " -D " + GATKDataLocation + "dbsnp_129_b36.rod" +
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
        new VCITTest("-L 1:1-10000 --printPerLocus", "e4ee2eaa3114888e918a1c82df7a027a");
        new VCITTest("-L 1:1-10000 --printPerLocus --onlyContextsOfType SNP", "2097e32988d603d3b353b50218c86d3b");
        new VCITTest("-L 1:1-10000 --printPerLocus --onlyContextsOfType INDEL", "033bd952fca048fe1a4f6422b57ab2ed");
        new VCITTest("-L 1:1-10000 --printPerLocus --onlyContextsOfType MIXED", "e5a00766f8c1ff9cf92310bafdec3126");
        new VCITTest("-L 1:1-10000 --printPerLocus --onlyContextsOfType NO_VARIATION", "39335acdb34c8a2af433dc50d619bcbc");

        // TODO : Eric, these are bad because the conversion fails
        //new VCITTest("-L 1:1-10000 --printPerLocus --takeFirstOnly", "5b5635e4877d82e8a27d70dac24bda2f");
        //new VCITTest("-L 1:1-10000 --printPerLocus --onlyContextsOfType INDEL --onlyContextsStartinAtCurrentPosition", "5e40980c02797f90821317874426a87a");
        //new VCITTest("-L 1:1-10000 --printPerLocus --onlyContextsStartinAtCurrentPosition", "ceced3f270b4fe407ee83bc9028becde");
        //new VCITTest("-L 1:1-10000 --printPerLocus --takeFirstOnly --onlyContextsStartinAtCurrentPosition", "9a9b9e283553c28bf58de1cafa38fe92");

        return VCITTest.getTests(VCITTest.class);
    }

    @Test(dataProvider = "VCITTestData")
    public void testConversionSelection(VCITTest test) {
	String extraArgs = test.args;
	String md5 = test.md5;

	WalkerTestSpec spec = new WalkerTestSpec( root + " " + extraArgs + " -o %s",
						  1, // just one output file
						  Arrays.asList(md5));
	executeTest("testDbSNPAndVCFConversions", spec);
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
                Arrays.asList("529f936aa6c303658b23caf4e527782f"));
         executeTest("testLargeScaleConversion", spec);
    }
}
