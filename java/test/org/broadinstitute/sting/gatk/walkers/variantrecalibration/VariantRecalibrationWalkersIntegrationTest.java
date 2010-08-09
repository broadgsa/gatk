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
        e.put( validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "05d1692624a28cd9446feac8fd2408ab" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String vcf = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -T GenerateVariantClusters" +
                            " -B input,VCF," + vcf +
                            " -L 1:1-100,000,000" +
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
        e.put( validationDataLocation + "yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.vcf", "df4e6f55c714b68b13e47b61d4bd0cd5" );
        
        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String vcf = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFiles.get(vcf);
            System.out.printf("PARAMS FOR %s is %s%n", vcf, paramsFile);
            if ( paramsFile != null ) {
                File file = createTempFile("cluster",".vcf");
                WalkerTestSpec spec = new WalkerTestSpec(
                        "-R " + b36KGReference +
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
