package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public class VariantRecalibrationWalkersIntegrationTest extends WalkerTest {
    private static class VRTest {
        String inVCF;
        String tranchesMD5;
        String recalMD5;
        String cutVCFMD5;
        public VRTest(String inVCF, String tranchesMD5, String recalMD5, String cutVCFMD5) {
            this.inVCF = inVCF;
            this.tranchesMD5 = tranchesMD5;
            this.recalMD5 = recalMD5;
            this.cutVCFMD5 = cutVCFMD5;
        }

        @Override
        public String toString() {
            return "VRTest{inVCF='" + inVCF +"'}";
        }
    }

    VRTest lowPass = new VRTest(validationDataLocation + "phase1.projectConsensus.chr20.raw.snps.vcf",
            "4d08c8eee61dd1bdea8c5765f34e41f0",  // tranches
            "ce396fe4045e020b61471f6737dff36e",  // recal file
            "4f59bd61be900b25c6ecedaa68b9c8de"); // cut VCF

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
                        " -U LENIENT_VCF_PROCESSING -o %s" +
                        " -tranchesFile " + getMd5DB().getMD5FilePath(params.tranchesMD5, null) +
                        " -recalFile " + getMd5DB().getMD5FilePath(params.recalMD5, null),
                Arrays.asList(params.cutVCFMD5));
        spec.disableShadowBCF(); // TODO -- enable when we support symbolic alleles
        executeTest("testApplyRecalibration-"+params.inVCF, spec);
    }

    VRTest bcfTest = new VRTest(privateTestDir + "vqsr.bcf_test.snps.unfiltered.bcf",
            "6a1eef4d02857dbb117a15420b5c0ce9",  // tranches
            "238366af66b05b6d21749e799c25353d",  // recal file
            "3928d6bc5007becf52312ade70f14c42"); // cut VCF

    @DataProvider(name = "VRBCFTest")
    public Object[][] createVRBCFTest() {
        return new Object[][]{ {bcfTest} };
        //return new Object[][]{ {yriTrio}, {lowPass} }; // Add hg19 chr20 trio calls here
    }

    @Test(dataProvider = "VRBCFTest")
    public void testVariantRecalibratorWithBCF(VRTest params) {
        //System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -resource:known=true,prior=10.0 " + GATKDataLocation + "dbsnp_132_b37.leftAligned.vcf" +
                        " -resource:truth=true,training=true,prior=15.0 " + comparisonDataLocation + "Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf" +
                        " -resource:training=true,truth=true,prior=12.0 " + comparisonDataLocation + "Validated/Omni2.5_chip/Omni25_sites_1525_samples.b37.vcf" +
                        " -T VariantRecalibrator" +
                        " -input " + params.inVCF +
                        " -L 20:10,000,000-20,000,000" +
                        " --no_cmdline_in_header" +
                        " -an AC " + // integer value
                        " -an QD -an ReadPosRankSum -an FS -an InbreedingCoeff " + // floats value
                        " -mG 2 "+
                        " -recalFile %s" +
                        " -tranchesFile %s",
                2,
                Arrays.asList("bcf", "txt"),
                Arrays.asList(params.recalMD5, params.tranchesMD5));
        executeTest("testVariantRecalibrator-"+params.inVCF, spec).getFirst();
    }

    @Test(dataProvider = "VRBCFTest", dependsOnMethods="testVariantRecalibratorWithBCF")
    public void testApplyRecalibrationWithBCF(VRTest params) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " -T ApplyRecalibration" +
                        " -L 20:10,000,000-20,000,000" +
                        " --no_cmdline_in_header" +
                        " -input " + params.inVCF +
                        " -U LENIENT_VCF_PROCESSING -o %s" +
                        " -tranchesFile " + getMd5DB().getMD5FilePath(params.tranchesMD5, null) +
                        " -recalFile " + getMd5DB().getMD5FilePath(params.recalMD5, null),
                Arrays.asList(params.cutVCFMD5));
        spec.disableShadowBCF();
        executeTest("testApplyRecalibration-"+params.inVCF, spec);
    }


    VRTest indelUnfiltered = new VRTest(
            validationDataLocation + "combined.phase1.chr20.raw.indels.unfiltered.sites.vcf", // all FILTERs as .
            "b7589cd098dc153ec64c02dcff2838e4",  // tranches
            "a04a9001f62eff43d363f4d63769f3ee",  // recal file
            "b2c6827be592c24a4692b1753edc7d23"); // cut VCF

    VRTest indelFiltered = new VRTest(
            validationDataLocation + "combined.phase1.chr20.raw.indels.filtered.sites.vcf", // all FILTERs as PASS
            "b7589cd098dc153ec64c02dcff2838e4",  // tranches
            "a04a9001f62eff43d363f4d63769f3ee",  // recal file
            "5d483fe1ba2ef36ee9e6c14cbd654706"); // cut VCF

    @DataProvider(name = "VRIndelTest")
    public Object[][] createTestVariantRecalibratorIndel() {
        return new Object[][]{ {indelUnfiltered}, {indelFiltered} };
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
                        " -U LENIENT_VCF_PROCESSING --no_cmdline_in_header" +
                        " -input " + params.inVCF +
                        " -o %s" +
                        " -tranchesFile " + getMd5DB().getMD5FilePath(params.tranchesMD5, null) +
                        " -recalFile " + getMd5DB().getMD5FilePath(params.recalMD5, null),
                Arrays.asList(params.cutVCFMD5));
        spec.disableShadowBCF(); // has to be disabled because the input VCF is missing LowQual annotation
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
                        " -input " + privateTestDir + "VQSR.mixedTest.input" +
                        " -o %s" +
                        " -tranchesFile " + privateTestDir + "VQSR.mixedTest.tranches" +
                        " -recalFile " + privateTestDir + "VQSR.mixedTest.recal",
                Arrays.asList("018b3a5cc7cf0cb5468c6a0c80ccaa8b"));
        executeTest("testApplyRecalibrationSnpAndIndelTogether", spec);
    }
}

