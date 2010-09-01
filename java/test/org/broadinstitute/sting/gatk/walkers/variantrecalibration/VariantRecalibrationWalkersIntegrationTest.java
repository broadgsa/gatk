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
        e.put( validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "cb3e8d9072d478243c3f9d0ee09fef3b" );
        e.put( validationDataLocation + "lowpass.N3.chr1.raw.0.vcf", "0d4fee916a886ca98c286d3b7fed0ff6" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String vcf = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " --DBSNP " + GATKDataLocation + "dbsnp_129_b36.rod" +
                            " -B:hapmap,VCF " + validationDataLocation + "CEU_hapmap_nogt_23.vcf" +
                            " -weightDBSNP 1.0 -weightHapMap 1.0" +
                            " -T GenerateVariantClusters" +
                            " -B:input,VCF " + vcf +
                            " -L 1:1-100,000,000" +
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
                Arrays.asList("ad0868adddb2e837b4fc08f140e4d9c3", "937080353c7e03e11f8a70fc0004bf76","5b89fa5a4edf0080d64230d4103d2b8d")); // Each test checks the md5 of three output files
        e.put( validationDataLocation + "lowpass.N3.chr1.raw.0.vcf",
                Arrays.asList("054d3228acfd6b02d24bfcf2fbd280a0", "f7c5c6cff9dd5280b25e24e0591e4cb0","1de1473db5720b882edf1381fa3dd039")); // Each test checks the md5 of three output files

        for ( Map.Entry<String, List<String>> entry : e.entrySet() ) {
            String vcf = entry.getKey();
            List<String> md5s = entry.getValue();
            String clusterFile = clusterFiles.get(vcf);
            System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
            if ( clusterFile != null ) {
                WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                        "-R " + b36KGReference +
                                " --DBSNP " + GATKDataLocation + "dbsnp_129_b36.rod" +
                                " -B:hapmap,VCF " + validationDataLocation + "CEU_hapmap_nogt_23.vcf" +
                                " -T VariantRecalibrator" +
                                " -B:input,VCF " + vcf +
                                " -L 1:40,000,000-100,000,000" +
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
        e.put( validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "841deca1cb30a5a081cfbaa6f663e22a" );
        e.put( validationDataLocation + "lowpass.N3.chr1.raw.0.vcf", "6234259e75fd7c1db31566dde1c55a82" );

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
                                " -L 1:40,000,000-100,000,000" +
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
