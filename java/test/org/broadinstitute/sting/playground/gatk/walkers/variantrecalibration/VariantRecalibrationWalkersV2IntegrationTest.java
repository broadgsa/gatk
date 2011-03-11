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
            "f65c27ee40053adc72dd0bfbb628e4d7",  // cluster file
            "dce581b880ffb6ea39cbada1ecc95915",  // tranches
            "c3e8a2f43656eab7d847dbf850f844a6",  // recalVCF
            "50f752a72643db9ad0aa94b3fc4e23d6"); // cut VCF

    VRTest lowPass = new VRTest("lowpass.N3.chr1.raw.vcf",
            "bda8f17cfc19d23e7e51f99e547f4b3d",  // cluster file
            "66edae83c50f4e8601fef7fafba774af",  // tranches
            "0123537e373657386068a534c0f5c91b",  // recalVCF
            "2172368e8585841e5ad96c95d0827c4b"); // cut VCF

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
                        " -B:hapmap,VCF " + comparisonDataLocation + "Validated/HapMap/3.2/sites_r27_nr.b36_fwd.vcf" +
                        " -weightDBSNP 1.0 -weightHapMap 1.0" +
                        " -T GenerateVariantClustersV2" +
                        " -B:input,VCF " + params.inVCF +
                        " -L 1:50,000,000-200,000,000" +
                        " -qual 50.0" +
                        " --ignore_filter GATK_STANDARD" +
                        " -an QD -an MQ -an SB" +
                        " -clusterFile %s",
                Arrays.asList(params.clusterMD5));
        executeTest("testGenerateVariantClusters-"+params.inVCF, spec).getFirst();
    }

    @Test(dataProvider = "VRTest",dependsOnMethods="testGenerateVariantClusters")
    public void testVariantRecalibrator(VRTest params) {
        //System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -NO_HEADER" +
                        " --DBSNP " + GATKDataLocation + "dbsnp_129_b36.rod" +
                        " -B:hapmap,VCF " + comparisonDataLocation + "Validated/HapMap/3.2/sites_r27_nr.b36_fwd.vcf" +
                        " -B:truthHapmap,VCF " + comparisonDataLocation + "Validated/HapMap/3.2/sites_r27_nr.b36_fwd.vcf" +
                        " -T VariantRecalibratorV2" +
                        " -B:input,VCF " + params.inVCF +
                        " -L 1:20,000,000-100,000,000" +
                        " --ignore_filter GATK_STANDARD" +
                        " --ignore_filter HARD_TO_VALIDATE" +
                        " -clusterFile " + getFileForMD5(params.clusterMD5) +
                        " -sm TRUTH_SENSITIVITY" +
                        " -o %s" +
                        " -tranchesFile %s",
                Arrays.asList(params.recalVCFMD5, params.tranchesMD5));
        executeTest("testVariantRecalibrator-"+params.inVCF, spec).getFirst();
    }

    @Test(dataProvider = "VRTest",dependsOnMethods="testVariantRecalibrator")
    public void testApplyVariantCuts(VRTest params) {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -NO_HEADER" +
                        " --DBSNP " + GATKDataLocation + "dbsnp_129_b36.rod" +
                        " -T ApplyVariantCutsV2" +
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
                        " -T GenerateVariantClustersV2" +
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

