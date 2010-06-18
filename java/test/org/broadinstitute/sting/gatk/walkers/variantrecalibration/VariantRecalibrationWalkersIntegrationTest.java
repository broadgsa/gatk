package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.*;
import java.io.File;

public class VariantRecalibrationWalkersIntegrationTest extends WalkerTest {
    static HashMap<String, String> paramsFiles = new HashMap<String, String>();
    
    @Test
    public void testGenerateVariantClusters() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "9b7517ff1fd0fc23a2596acc82e2ed96" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String vcf = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -T GenerateVariantClusters" +
                            " -B input,VCF," + vcf +
                            " -L 1:1-100,000,000" +
                            " -nG 6" +
                            " --ignore_filter GATK_STANDARD" +
                            " -an QD -an HRun -an SB" +
                            " -clusterFile %s",
                    1, // just one output file
                    Arrays.asList(md5));
            List<File> result = executeTest("testGenerateVariantClusters", spec).getFirst();
            paramsFiles.put(vcf, result.get(0).getAbsolutePath());
        }
    }

    @Test
    public void testVariantRecalibrator() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "e0e3d959929aa3940a81c9926d1406e2" );
        
        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String vcf = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFiles.get(vcf);
            System.out.printf("PARAMS FOR %s is %s%n", vcf, paramsFile);
            if ( paramsFile != null ) {
                File file = createTempFile("cluster",".vcf");
                WalkerTestSpec spec = new WalkerTestSpec(
                        "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                                " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                                " -T VariantRecalibrator" +
                                " -B input,VCF," + vcf +
                                " -L 1:40,000,000-100,000,000" +
                                " --ignore_filter GATK_STANDARD" +
                                " -output " + file.getAbsolutePath().substring(0,file.getAbsolutePath().length()-4) +
                                " -clusterFile " + paramsFile,
                        0,
                        new ArrayList<String>(0));

                spec.addAuxFile(md5, file);

                executeTest("testVariantRecalibrator", spec);
            }
        }
    }
}
