package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
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
                            " -NO_HEADER" +
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
                Arrays.asList("ca80c95e47d22a79eb0700c039457aca", "ab6b91a3f97b682817d2cc068c34c317","a4802dbf883138ae1152109be990fb9d")); // Each test checks the md5 of three output files
        e.put( validationDataLocation + "lowpass.N3.chr1.raw.vcf",
                Arrays.asList("a40d3b1dd7bfbb66a52145600f87d744", "e2c89c74debc1c202c56060b77575dff","8e8b55b521aaba1c44169a19e3ff2355")); // Each test checks the md5 of three output files

        for ( Map.Entry<String, List<String>> entry : e.entrySet() ) {
            String vcf = entry.getKey();
            List<String> md5s = entry.getValue();
            String clusterFile = clusterFiles.get(vcf);
            System.out.printf("PARAMS FOR %s is %s%n", vcf, clusterFile);
            if ( clusterFile != null ) {
                WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                        "-R " + b36KGReference +
                                " -NO_HEADER" +
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
                        3, // three output files
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
        e.put( validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "ece0c15c34926fc585e12503f6ce6271" );
        e.put( validationDataLocation + "lowpass.N3.chr1.raw.vcf", "994b329a35d01e9564f5581cf3d9feac" );

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
                                " -NO_HEADER" +
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


    @Test
    public void testFailWithBadAnnotation() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "lowpass.N3.chr1.raw.vcf", "" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String vcf = entry.getKey();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -NO_HEADER" +
                            " --DBSNP " + GATKDataLocation + "dbsnp_129_b36.rod" +
                            " -B:hapmap,VCF " + comparisonDataLocation + "Validated/HapMap/3.2/genotypes_r27_nr.b36_fwd.vcf" +
                            " -weightDBSNP 0.2 -weightHapMap 1.0" +
                            " -T GenerateVariantClusters" +
                            " -B:input,VCF " + vcf +
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
}

