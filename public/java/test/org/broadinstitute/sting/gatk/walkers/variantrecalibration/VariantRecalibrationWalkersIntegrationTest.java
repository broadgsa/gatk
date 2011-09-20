package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.broadinstitute.sting.MD5DB;
import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.util.*;

public class VariantRecalibrationWalkersIntegrationTest extends WalkerTest {
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

    VRTest lowPass = new VRTest("phase1.projectConsensus.chr20.raw.snps.vcf",
            "0ddd1e0e483d2eaf56004615cea23ec7",  // tranches
            "58780f63182e139fdbe17f6c18b5b774",  // recal file
            "f67d844b6252a55452cf4167b77530b1"); // cut VCF

    @DataProvider(name = "VRTest")
    public Object[][] createData1() {
        return new Object[][]{ {lowPass} };
        //return new Object[][]{ {yriTrio}, {lowPass} }; // Add hg19 chr20 trio calls here
    }

    @Test(dataProvider = "VRTest")
    public void testVariantRecalibrator(VRTest params) {
        //System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -resource:known=true,prior=10.0 " + GATKDataLocation + "dbsnp_132_b37.leftAligned.vcf" +
                        " -resource:truth=true,training=true,prior=15.0 " + comparisonDataLocation + "Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf" +
                        " -resource:training=true,truth=true,prior=12.0 " + comparisonDataLocation + "Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf" +
                        " -T VariantRecalibrator" +
                        " -input " + params.inVCF +
                        " -L 20:1,000,000-40,000,000" +
                        " -an QD -an HaplotypeScore -an HRun" +
                        " -percentBad 0.07" +
                        " --minNumBadVariants 0" +
                        " --trustAllPolymorphic" + // for speed
                        " -recalFile %s" +
                        " -tranchesFile %s",
                Arrays.asList(params.recalMD5, params.tranchesMD5));
        executeTest("testVariantRecalibrator-"+params.inVCF, spec).getFirst();
    }

    @Test(dataProvider = "VRTest",dependsOnMethods="testVariantRecalibrator")
    public void testApplyRecalibration(VRTest params) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -T ApplyRecalibration" +
                        " -L 20:12,000,000-30,000,000" +
                        " -NO_HEADER" +
                        " -input " + params.inVCF +
                        " -o %s" +
                        " -tranchesFile " + MD5DB.getMD5FilePath(params.tranchesMD5, null) +
                        " -recalFile " + MD5DB.getMD5FilePath(params.recalMD5, null),
                Arrays.asList(params.cutVCFMD5));
        executeTest("testApplyRecalibration-"+params.inVCF, spec);
    }

    VRTest indel = new VRTest("combined.phase1.chr20.raw.indels.sites.vcf",
            "6d7ee4cb651c8b666e4a4523363caaff",  // tranches
            "4759b111a5aa53975d46e0f22c7983bf",  // recal file
            "5d7e07d8813db96ba3f3dfe4737f83d1"); // cut VCF

    @DataProvider(name = "VRIndelTest")
    public Object[][] createData2() {
        return new Object[][]{ {indel} };
    }

    @Test(dataProvider = "VRIndelTest")
    public void testVariantRecalibratorIndel(VRTest params) {
        //System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -resource:known=true,prior=10.0 " + GATKDataLocation + "dbsnp_132_b37.leftAligned.vcf" +
                        " -resource:training=true,truth=true,prior=15.0 " + comparisonDataLocation + "Validated/Mills_Devine_Indels_2011/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.vcf" +
                        " -T VariantRecalibrator" +
                        " -input " + params.inVCF +
                        " -L 20:1,000,000-40,000,000" +
                        " -an QD -an ReadPosRankSum -an HaplotypeScore" +
                        " -percentBad 0.08" +
                        " -mode INDEL -mG 3" +
                        " --minNumBadVariants 0" +
                        " --trustAllPolymorphic" + // for speed
                        " -recalFile %s" +
                        " -tranchesFile %s",
                Arrays.asList(params.recalMD5, params.tranchesMD5));
        executeTest("testVariantRecalibratorIndel-"+params.inVCF, spec).getFirst();
    }

    @Test(dataProvider = "VRIndelTest",dependsOnMethods="testVariantRecalibratorIndel")
    public void testApplyRecalibrationIndel(VRTest params) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -T ApplyRecalibration" +
                        " -L 20:12,000,000-30,000,000" +
                        " -mode INDEL" +
                        " -NO_HEADER" +
                        " -input " + params.inVCF +
                        " -o %s" +
                        " -tranchesFile " + MD5DB.getMD5FilePath(params.tranchesMD5, null) +
                        " -recalFile " + MD5DB.getMD5FilePath(params.recalMD5, null),
                Arrays.asList(params.cutVCFMD5));
        executeTest("testApplyRecalibrationIndel-"+params.inVCF, spec);
    }
}

