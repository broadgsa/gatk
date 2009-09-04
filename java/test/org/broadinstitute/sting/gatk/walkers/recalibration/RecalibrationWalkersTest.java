package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;
import java.util.List;
import java.io.File;

public class RecalibrationWalkersTest extends WalkerTest {
    static HashMap<String, String> paramsFiles = new HashMap<String, String>();

    @Test
    public void testCountCovariates1() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.SLX.SRP000031.2009_06.chr1.10_20mb.bam", "47664c48992f593258932583576b47e4" );
        //e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12762.SOLID.SRP000031.2009_07.chr1.10_20mb.bam", "" );
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "b58185dc5fbdd88ca9539d940dff6c1a" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R /broad/1KG/reference/human_b36_both.fasta" +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -L 1:10,000,000-11,000,000" +
                            " --params %s",
                    1, // just one output file
                    Arrays.asList(md5));
            List<File> result = executeTest("testCountCovariates1", spec).getFirst();
            paramsFiles.put(bam, result.get(0).getAbsolutePath());
        }
    }

    @Test
    public void testTableRecalibrator1() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.SLX.SRP000031.2009_06.chr1.10_20mb.bam", "c98525aca6493179f084159df0264782" );
        //e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12762.SOLID.SRP000031.2009_07.chr1.10_20mb.bam", "" );
        //e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "28c5026359a952a5c652b1ccc2acac04" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFiles.get(bam);
            System.out.printf("PARAMS FOR %s is %s%n", bam, paramsFile);

            WalkerTestSpec spec = new WalkerTestSpec(
                    "-R /broad/1KG/reference/human_b36_both.fasta" +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -T TableRecalibration" +
                            " -I " + bam +
                            " -L 1:10,000,000-20,000,000" +
                            " --outputBam %s" +
                            " --params " + paramsFile,
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testTableRecalibrator1", spec);
        }
    }
}