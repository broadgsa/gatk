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
        e.put( validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "c8a0eaeed9a4f8c12b90d89c65ad3405" );
        e.put( validationDataLocation + "lowpass.N3.chr1.raw.vcf", "ead6836dcc9fde2dd26f42317395e92d" );

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
                Arrays.asList("274b58b7e45619411b061610d2ae0b3f", "e96c86d63be4401414dbadd43181e433","acd93e4747f5abb6ab81755a229168b5")); // Each test checks the md5 of three output files
        e.put( validationDataLocation + "lowpass.N3.chr1.raw.vcf",
                Arrays.asList("80176c57b9d9bbacee36e8dec9ec2c93", "f12879c5d4aa214c50d50f2d0fd6a60b","59f0f0d2bdb9f2f3c112a71f38dc834f")); // Each test checks the md5 of three output files

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
        e.put( validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "52538fbc9113271e1a0bf0ae3c904c93" );
        e.put( validationDataLocation + "lowpass.N3.chr1.raw.vcf", "b70686bf091d0306fd23daef729d05fb" );

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
