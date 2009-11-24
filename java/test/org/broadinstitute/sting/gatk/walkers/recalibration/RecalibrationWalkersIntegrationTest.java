package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;
import java.util.List;
import java.io.File;

public class RecalibrationWalkersIntegrationTest extends WalkerTest {
    static HashMap<String, String> paramsFiles = new HashMap<String, String>();

    @Test
    public void testCountCovariates1() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12892.SLX.SRP000031.2009_06.selected.bam", "5bb3cb7445ec223cd450b05430cb00cd" );
        //e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12762.SOLID.SRP000031.2009_07.chr1.10_20mb.bam", "dc1a1b6e99f4a47525cc1dce7b6eb1dc" ); // This file doesn't have read groups?
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "7842b33245a53b13c7b6cd4e4616ac11");
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "cbe17568c5f03516a88ff0c5ce33a87b" );
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "8806fb7f0e669d1abab61dc4b3d10d8f" );
        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R /broad/1KG/reference/human_b36_both.fasta" +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            ( bam.equals( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" )
                                    ? " -L 1:10,800,000-10,810,000" : " -L 1:10,000,000-10,200,000" ) +
                            " -cov ReadGroupCovariate" +
                            " -cov QualityScoreCovariate" +
                            " -cov CycleCovariate" +
                            " -cov DinucCovariate" +
                            " --sorted_output" +
                            " -recalFile %s",
                    1, // just one output file
                    Arrays.asList(md5));
            List<File> result = executeTest("testCountCovariates1", spec).getFirst();
            paramsFiles.put(bam, result.get(0).getAbsolutePath());
        }
    }

    @Test
    public void testTableRecalibrator1() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12892.SLX.SRP000031.2009_06.selected.bam", "d7f3e0db5ed9fefc917144a0f937d50d" );
        //e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12762.SOLID.SRP000031.2009_07.chr1.10_20mb.bam", "7be0b7a624d22187e712131f12aae546" );
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "94f2c602fef9800e270326be9faab3ad");
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "eea2e02ffe71c9a1f318d2e9cd31a103" );
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "8126b9494837d504ef1277a799267c15" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFiles.get(bam);
            System.out.printf("PARAMS FOR %s is %s%n", bam, paramsFile);
            if ( paramsFile != null ) {
                WalkerTestSpec spec = new WalkerTestSpec(
                        "-R /broad/1KG/reference/human_b36_both.fasta" +
                                " -T TableRecalibration" +
                                " -I " + bam +
                                ( bam.equals( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" )
                                    ? " -L 1:10,800,000-10,810,000" : " -L 1:10,100,000-10,300,000" ) +
                                " -outputBam %s" +
                                " --no_pg_tag" +
                                " -recalFile " + paramsFile,
                        1, // just one output file
                        Arrays.asList(md5));
                executeTest("testTableRecalibrator1", spec);
            }
        }
    }
}
