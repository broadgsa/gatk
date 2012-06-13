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
            "6e1f98bb819ccf03e17a2288742160d3",  // recal file
            "1050c387d170639f8cec221e5dddd626"); // cut VCF

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
                        " --no_cmdline_in_header" +
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
                        " --no_cmdline_in_header" +
                        " -input " + params.inVCF +
                        " -o %s" +
                        " -tranchesFile " + getMd5DB().getMD5FilePath(params.tranchesMD5, null) +
                        " -recalFile " + getMd5DB().getMD5FilePath(params.recalMD5, null),
                Arrays.asList(params.cutVCFMD5));
        executeTest("testApplyRecalibration-"+params.inVCF, spec);
    }

    VRTest indel = new VRTest("combined.phase1.chr20.raw.indels.sites.vcf",
            "da4458d05f6396f5c4ab96f274e5ccdc",  // tranches
            "8e2417336fa62e6c4d9f61b6deebdd82",  // recal file
            "bf0e8ed5e250d52f0545074c61217d16"); // cut VCF

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
                        " --no_cmdline_in_header" +
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
                        " --no_cmdline_in_header" +
                        " -input " + params.inVCF +
                        " -o %s" +
                        " -tranchesFile " + getMd5DB().getMD5FilePath(params.tranchesMD5, null) +
                        " -recalFile " + getMd5DB().getMD5FilePath(params.recalMD5, null),
                Arrays.asList(params.cutVCFMD5));
        executeTest("testApplyRecalibrationIndel-"+params.inVCF, spec);
    }

    @Test
    public void testApplyRecalibrationSnpAndIndelTogether() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -T ApplyRecalibration" +
                        " -L 20:1000100-1000500" +
                        " -mode BOTH" +
                        " --no_cmdline_in_header" +
                        " -input " + testDir + "VQSR.mixedTest.input" +
                        " -o %s" +
                        " -tranchesFile " + testDir + "VQSR.mixedTest.tranches" +
                        " -recalFile " + testDir + "VQSR.mixedTest.recal",
                Arrays.asList("1370d7701a6231633d43a8062b7aff7f"));
        executeTest("testApplyRecalibrationSnpAndIndelTogether", spec);
    }
}

