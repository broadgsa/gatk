package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.util.*;
import java.io.File;

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
            "aa211561e6e9b1e323dec4eeaa318088",  // tranches
            "226adf939e90ad20599108e77ad25dae",  // recal file
            "345a159fd54f5c38500bb20f2de13737"); // cut VCF

    VRTest lowPass = new VRTest("lowpass.N3.chr1.raw.vcf",
            "f4f2057a8c010206f0f56deff0602452",  // tranches
            "dd36d252a6dc6e3207addddae731dd88",  // recal file
            "f1ffd3bb1da3863ccb298a3373d6590a"); // cut VCF

    @DataProvider(name = "VRTest")
    public Object[][] createData1() {
        return new Object[][]{ {yriTrio}, {lowPass} };
    }

    @Test(dataProvider = "VRTest")
    public void testVariantRecalibrator(VRTest params) {
        //System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --DBSNP " + GATKDataLocation + "dbsnp_129_b36.rod" +
                        " -B:hapmap,VCF " + comparisonDataLocation + "Validated/HapMap/3.2/sites_r27_nr.b36_fwd.vcf" +
                        " -B:omni,VCF " + comparisonDataLocation + "Validated/Omni2.5_chip/1212samples.b36.vcf" +
                        " -T ContrastiveRecalibrator" +
                        " -B:input,VCF " + params.inVCF +
                        " -L 1:50,000,000-120,000,000" +
                        " -an QD -an MQ -an SB -an AF" +
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
                        " -B:input,VCF " + params.inVCF +
                        " -o %s" +
                        " -tranchesFile " + getFileForMD5(params.tranchesMD5) +
                        " -recalFile " + getFileForMD5(params.recalMD5),
                Arrays.asList(params.cutVCFMD5));
        executeTest("testApplyRecalibration-"+params.inVCF, spec);
    }
}

