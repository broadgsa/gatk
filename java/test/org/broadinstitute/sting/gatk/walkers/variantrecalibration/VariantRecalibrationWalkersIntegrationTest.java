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
            "4eeffa7a1965ce0c25c5edd0bae76290",  // in vcf
            "7407987a0148284ed910e1858116dd8d",  // tranches
            "15ab55be5b2f62627aea8546a4728d77",  // recalVCF
            "9435f1aed7313fbfff540a4d6d19d0c4"); // cut VCF

    VRTest lowPass = new VRTest("lowpass.N3.chr1.raw.vcf",
            "8937a3ae7f176dacf47b8ee6c0023416",  // in vcf
            "2896657b5c30bfd8e82e62e58d94ef4e",  // tranches
            "a5fe2ee50144ef61121c42daf430381c",  // recalVCF
            "9a35b69bed93894306c87bc9a0bcc116"); // cut VCF

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

    @Test(dataProvider = "VRTest",dependsOnMethods="testGenerateVariantClusters")
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

    @Test(dataProvider = "VRTest",dependsOnMethods="testVariantRecalibrator")
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

