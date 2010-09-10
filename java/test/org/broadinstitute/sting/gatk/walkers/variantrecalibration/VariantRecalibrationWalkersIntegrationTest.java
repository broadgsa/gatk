package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.*;
import java.io.File;

public class VariantRecalibrationWalkersIntegrationTest extends WalkerTest {
    static HashMap<String, String> clusterFiles = new HashMap<String, String>();
    static HashMap<String, String> tranchesFiles = new HashMap<String, String>();
    static HashMap<String, String> inputVCFFiles = new HashMap<String, String>();

    @Test
    public void testGenerateVariantClusters() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "ab2629d67e378fd3aceb8318f0fbfe04" );
        e.put( validationDataLocation + "lowpass.N3.chr1.raw.vcf", "725489156426e4ddd8d623ab3d4b1023" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String vcf = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " --DBSNP " + GATKDataLocation + "dbsnp_129_b36.rod" +
                            " -B:hapmap,VCF " + comparisonDataLocation + "Validated/HapMap/3.2/genotypes_r27_nr.b36_fwd.vcf" +
                            " -weightDBSNP 0.2 -weightHapMap 1.0" +
                            " -T GenerateVariantClusters" +
                            " -B:input,VCF " + vcf +
                            " -L 1:50,000,000-200,000,000" +
                            " -qual 50.0" +
                            " --ignore_filter GATK_STANDARD" +
                            " -an QD -an HRun -an SB" +
                            " -clusterFile %s",
                    1, // just one output file
                    Arrays.asList(md5));
            List<File> result = executeTest("testGenerateVariantClusters", spec).getFirst();
            clusterFiles.put(vcf, result.get(0).getAbsolutePath());
        }
    }

    @Test
    public void testVariantRecalibrator() {
        HashMap<String, List<String>> e = new HashMap<String, List<String>>();
        e.put( validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf",
                Arrays.asList("e7dc06c6d4d8c2cdd03bfdac483aa92b", "038c31c5bb46a4df89b8ee69ec740812","7d42bbdfb69fdfb18cbda13a63d92602")); // Each test checks the md5 of three output files
        e.put( validationDataLocation + "lowpass.N3.chr1.raw.vcf",
                Arrays.asList("80243bd6731f6b0341fa1762c095a72e", "661360e85392af9c97e386399871854a","371e5a70a4006420737c5ab259e0e23e")); // Each test checks the md5 of three output files

        for ( Map.Entry<String, List<String>> entry : e.entrySet() ) {
            String vcf = entry.getKey();
            List<String> md5s = entry.getValue();
            String clusterFile = clusterFiles.get(vcf);
            System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
            if ( clusterFile != null ) {
                WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                        "-R " + b36KGReference +
                                " --DBSNP " + GATKDataLocation + "dbsnp_129_b36.rod" +
                                " -B:hapmap,VCF " + comparisonDataLocation + "Validated/HapMap/3.2/genotypes_r27_nr.b36_fwd.vcf" +
                                " -T VariantRecalibrator" +
                                " -B:input,VCF " + vcf +
                                " -L 1:20,000,000-100,000,000" +
                                " --ignore_filter GATK_STANDARD" +
                                " --ignore_filter HARD_TO_VALIDATE" +
                                " -clusterFile " + clusterFile +
                                " -titv 2.07" +
                                " -qScale 20.0" + 
                                " -o %s" +
                                " -tranchesFile %s" +
                                " -reportDatFile %s",                                
                        3, // two output file
                        md5s);
                List<File> result = executeTest("testVariantRecalibrator", spec).getFirst();
                inputVCFFiles.put(vcf, result.get(0).getAbsolutePath());
                tranchesFiles.put(vcf, result.get(1).getAbsolutePath());
            }
        }


    }

    @Test
    public void testApplyVariantCuts() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "ee22087a813aadf7d9426cab5a6a2164" );
        e.put( validationDataLocation + "lowpass.N3.chr1.raw.vcf", "d9e00e6fe8269b3218880d2c084804c6" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String vcf = entry.getKey();
            String md5 = entry.getValue();
            String inputVCFFile = inputVCFFiles.get(vcf);
            String tranchesFile = tranchesFiles.get(vcf);
            System.out.printf("PARAMS FOR %s is %s%n", vcf, inputVCFFile);
            System.out.printf("PARAMS FOR %s is %s%n", vcf, tranchesFile);
            if ( inputVCFFile != null && tranchesFile != null ) {
                WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                        "-R " + b36KGReference +
                                " --DBSNP " + GATKDataLocation + "dbsnp_129_b36.rod" +
                                " -T ApplyVariantCuts" +
                                " -L 1:20,000,000-100,000,000" +
                                " -B:input,VCF " + inputVCFFile +
                                " -o %s" +
                                " -tranchesFile " + tranchesFile,
                        1, // just one output file
                        Arrays.asList(md5));
                executeTest("testApplyVariantCuts", spec);
            }
        }
    }
}
