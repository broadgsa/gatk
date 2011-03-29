package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.util.*;

public class VariantRecalibrationWalkersV2IntegrationTest extends WalkerTest {
    static HashMap<String, String> clusterFiles = new HashMap<String, String>();
    static HashMap<String, String> tranchesFiles = new HashMap<String, String>();
    static HashMap<String, String> inputVCFFiles = new HashMap<String, String>();

    private static class VRTest {
        String inVCF;
        String tranchesMD5;
        String recalMD5;
        String cutVCFMD5;
        public VRTest(String inVCF, String tranchesMD5, String recalMD5, String cutVCFMD5) {
            this.inVCF = validationDataLocation + inVCF;
            this.tranchesMD5 = tranchesMD5;
            this.recalMD5 = recalMD5;
            this.cutVCFMD5 = cutVCFMD5;
        }
    }

    VRTest yriTrio = new VRTest("yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf",
            "644e37fe6edaee13ddf9f3c4f0e6f500",  // tranches
            "e467bfd4ab74c7a7ff303cf810c514b3",  // recal file
            "05f66c959d1f8e78bcb4298092524921"); // cut VCF

    VRTest lowPass = new VRTest("lowpass.N3.chr1.raw.vcf",
            "e0a3845619d2138717df39ef1e9ee75f",  // tranches
            "c3fa524d0bdb40d6bd7aeda9e8dd882c",  // recal file
            "9f21a668c3fff13f382367851a11303c"); // cut VCF

    @DataProvider(name = "VRTest")
    public Object[][] createData1() {
        return new Object[][]{ {yriTrio}, {lowPass} };
    }

    @Test(dataProvider = "VRTest")
    public void testVariantRecalibrator(VRTest params) {
        //System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -B:dbsnp,DBSNP,known=true,training=false,truth=false,prior=10.0 " + GATKDataLocation + "dbsnp_129_b36.rod" +
                        " -B:hapmap,VCF,known=false,training=true,truth=true,prior=12.0 " + comparisonDataLocation + "Validated/HapMap/3.2/sites_r27_nr.b36_fwd.vcf" +
                        " -B:omni,VCF,known=false,training=true,truth=true,prior=15.0 " + comparisonDataLocation + "Validated/Omni2.5_chip/1212samples.b36.vcf" +
                        " -T ContrastiveRecalibrator" +
                        " -B:input,VCF " + params.inVCF +
                        " -L 1:50,000,000-120,000,000" +
                        " -an QD -an MQ -an SB -an AF" +
                        " --trustAllPolymorphic" + // for speed
                        " -recalFile %s" +
                        " -tranchesFile %s",
                Arrays.asList(params.recalMD5, params.tranchesMD5));
        executeTest("testVariantRecalibrator-"+params.inVCF, spec).getFirst();
    }

    @Test(dataProvider = "VRTest",dependsOnMethods="testVariantRecalibrator")
    public void testApplyRecalibration(VRTest params) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -T ApplyRecalibration" +
                        " -L 1:60,000,000-115,000,000" +
                        " -NO_HEADER" +
                        " -B:input,VCF " + params.inVCF +
                        " -o %s" +
                        " -tranchesFile " + getFileForMD5(params.tranchesMD5) +
                        " -recalFile " + getFileForMD5(params.recalMD5),
                Arrays.asList(params.cutVCFMD5));
        executeTest("testApplyRecalibration-"+params.inVCF, spec);
    }
}

