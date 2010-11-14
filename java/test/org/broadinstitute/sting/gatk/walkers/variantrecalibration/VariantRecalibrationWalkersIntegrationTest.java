package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.util.*;
import java.io.File;

public class VariantRecalibrationWalkersIntegrationTest extends WalkerTest {
    static HashMap<String, String> clusterFiles = new HashMap<String, String>();
    static HashMap<String, String> tranchesFiles = new HashMap<String, String>();
    static HashMap<String, String> inputVCFFiles = new HashMap<String, String>();

    private static class VRTest {
        String inVCF;
        String clusterMD5;
        String tranchesMD5;
        String recalVCFMD5;
        String cutVCFMD5;
        public VRTest(String inVCF, String clusterMD5, String tranchesMD5, String recalVCFMD5, String cutVCFMD5) {
            this.inVCF = validationDataLocation + inVCF;
            this.clusterMD5 = clusterMD5;
            this.tranchesMD5 = tranchesMD5;
            this.recalVCFMD5 = recalVCFMD5;
            this.cutVCFMD5 = cutVCFMD5;
        }
    }

    VRTest yriTrio = new VRTest("yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf",
            "ab2629d67e378fd3aceb8318f0fbfe04",  // in vcf
            "4dd95e9d8e5d21a6ab73de67d9663492",  // tranches
            "4e893672230fca625f70b0491f3b36cb",  // recalVCF
            "371b0e2796982485ae050e46892f6660"); // cut VCF

    VRTest lowPass = new VRTest("lowpass.N3.chr1.raw.vcf",
            "725489156426e4ddd8d623ab3d4b1023",  // in vcf
            "dfc07132f592a811d0d6c25a4cb67a09",  // tranches
            "d52e4f511c9c00f8c21dffea81c47103",  // recalVCF
            "39716b2f03b50e88855d7975dd1b9b3e"); // cut VCF

    @DataProvider(name = "VRTest")
    public Object[][] createData1() {
        return new Object[][]{ {yriTrio}, {lowPass} };
    }

    @Test(dataProvider = "VRTest", enabled = true)
    public void testGenerateVariantClusters(VRTest params) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -NO_HEADER" +
                        " --DBSNP " + GATKDataLocation + "dbsnp_129_b36.rod" +
                        " -B:hapmap,VCF " + comparisonDataLocation + "Validated/HapMap/3.2/genotypes_r27_nr.b36_fwd.vcf" +
                        " -weightDBSNP 0.2 -weightHapMap 1.0" +
                        " -T GenerateVariantClusters" +
                        " -B:input,VCF " + params.inVCF +
                        " -L 1:50,000,000-200,000,000" +
                        " -qual 50.0" +
                        " --ignore_filter GATK_STANDARD" +
                        " -an QD -an HRun -an SB" +
                        " -clusterFile %s",
                Arrays.asList(params.clusterMD5));
        executeTest("testGenerateVariantClusters-"+params.inVCF, spec).getFirst();
    }

    @Test(dataProvider = "VRTest", dependsOnMethods = {"testGenerateVariantClusters"}, enabled = true)
    public void testVariantRecalibrator(VRTest params) {
        //System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -NO_HEADER" +
                        " --DBSNP " + GATKDataLocation + "dbsnp_129_b36.rod" +
                        " -B:hapmap,VCF " + comparisonDataLocation + "Validated/HapMap/3.2/sites_r27_nr.b36_fwd.vcf" +
                        " -T VariantRecalibrator" +
                        " -B:input,VCF " + params.inVCF +
                        " -L 1:20,000,000-100,000,000" +
                        " --ignore_filter GATK_STANDARD" +
                        " --ignore_filter HARD_TO_VALIDATE" +
                        " -clusterFile " + getFileForMD5(params.clusterMD5) +
                        " -titv 2.07" +
                        " -o %s" +
                        " -tranchesFile %s",
                Arrays.asList(params.recalVCFMD5, params.tranchesMD5));
        executeTest("testVariantRecalibrator-"+params.inVCF, spec).getFirst();
    }

    @Test(dataProvider = "VRTest", dependsOnMethods = {"testVariantRecalibrator"}, enabled = true)
    public void testApplyVariantCuts(VRTest params) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -NO_HEADER" +
                        " --DBSNP " + GATKDataLocation + "dbsnp_129_b36.rod" +
                        " -T ApplyVariantCuts" +
                        " -L 1:20,000,000-100,000,000" +
                        " -B:input,VCF " + getFileForMD5(params.recalVCFMD5) +
                        " -o %s" +
                        " -tranchesFile " + getFileForMD5(params.tranchesMD5),
                Arrays.asList(params.cutVCFMD5));
        executeTest("testApplyVariantCuts-"+params.inVCF, spec);
    }


    @Test()
    public void testFailWithBadAnnotation() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -NO_HEADER" +
                        " --DBSNP " + GATKDataLocation + "dbsnp_129_b36.rod" +
                        " -B:hapmap,VCF " + comparisonDataLocation + "Validated/HapMap/3.2/sites_r27_nr.b36_fwd.vcf" +
                        " -weightDBSNP 0.2 -weightHapMap 1.0" +
                        " -T GenerateVariantClusters" +
                        " -B:input,VCF " + lowPass.inVCF +
                        " -L 1:50,000,000-200,000,000" +
                        " -qual 50.0" +
                        " --ignore_filter GATK_STANDARD" +
                        " -an QD -an HRun -an ThisAnnotationIsBAD" + // There is a bad annotation here
                        " -clusterFile %s",
                1, // just one output file
                UserException.MalformedFile.class);
        executeTest("testFailWithBadAnnotation", spec);
    }
}

