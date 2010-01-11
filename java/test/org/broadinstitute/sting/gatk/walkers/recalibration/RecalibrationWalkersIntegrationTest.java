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
    static HashMap<String, String> paramsFilesNoReadGroupTest = new HashMap<String, String>();

    @Test
    public void testCountCovariates1() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "e5b2d5a2f4283718dae678cbc84be847" );
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "89084b43b824f9e3c5e2afdfe0930542");
        e.put( validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "7d6428a76e07ed4b99351aa4df89634d" );
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "6ede6fc840c4e5070a58a919b48e7504" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            ( bam.equals( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" )
                                ? " -L 1:10,800,000-10,810,000" : " -L 1:10,000,000-10,200,000" ) +
                            " -cov ReadGroupCovariate" +
                            " -cov QualityScoreCovariate" +
                            " -cov CycleCovariate" +
                            " -cov DinucCovariate" +
                            " -cov TileCovariate" +
                            " --solid_recal_mode SET_Q_ZERO" +
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
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "6c59d291c37d053e0f188b762f3060a5" );
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "e06f1397b9c40f75e96cd3df76730ee0");
        e.put( validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "7ebdce416b72679e1cf88cc9886a5edc" );
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "48ddc93cae054f9423f3a7ed9f36540e" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFiles.get(bam);
            System.out.printf("PARAMS FOR %s is %s%n", bam, paramsFile);
            if ( paramsFile != null ) {
                WalkerTestSpec spec = new WalkerTestSpec(
                        "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                                " -T TableRecalibration" +
                                " -I " + bam +
                                ( bam.equals( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" )
                                    ? " -L 1:10,800,000-10,810,000" : " -L 1:10,100,000-10,300,000" ) +
                                " -outputBam %s" +
                                " --no_pg_tag" +
                                " --solid_recal_mode SET_Q_ZERO" +
                                " -recalFile " + paramsFile,
                        1, // just one output file
                        Arrays.asList(md5));
                executeTest("testTableRecalibrator1", spec);
            }
        }
    }

    @Test
    public void testTableRecalibratorMaxQ70() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "665711dfb81d67582b28faea24e26160" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFiles.get(bam);
            System.out.printf("PARAMS FOR %s is %s%n", bam, paramsFile);
            if ( paramsFile != null ) {
                WalkerTestSpec spec = new WalkerTestSpec(
                        "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                                " -T TableRecalibration" +
                                " -I " + bam +
                                ( bam.equals( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" )
                                    ? " -L 1:10,800,000-10,810,000" : " -L 1:10,100,000-10,300,000" ) +
                                " -outputBam %s" +
                                " --no_pg_tag" +
                                " -maxQ 70" +
                                " --solid_recal_mode SET_Q_ZERO" +
                                " -recalFile " + paramsFile,
                        1, // just one output file
                        Arrays.asList(md5));
                executeTest("testTableRecalibratorMaxQ70", spec);
            }
        }
    }

    @Test
    public void testCountCovariatesVCF() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "9b9d21ffb70f15ef2aebad21f3fc05cb");

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                            " -B dbsnp,VCF," + validationDataLocation + "vcfexample3.vcf" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -L 1:10,000,000-10,200,000" +
                            " -cov ReadGroupCovariate" +
                            " -cov QualityScoreCovariate" +
                            " -cov CycleCovariate" +
                            " -cov DinucCovariate" +
                            " --solid_recal_mode SET_Q_ZERO" +
                            " -recalFile %s",
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testCountCovariatesVCF", spec);
        }
    }

    @Test
    public void testCountCovariatesNoReadGroups() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12762.SOLID.SRP000031.2009_07.chr1.10_20mb.bam", "a86c64f649b847b7f81ac50a808d3d45" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -L 1:10,000,000-10,200,000" +
                            " -cov ReadGroupCovariate" +
                            " -cov QualityScoreCovariate" +
                            " -cov CycleCovariate" +
                            " -cov DinucCovariate" +
                            " --default_platform illumina" +
                            " --solid_recal_mode SET_Q_ZERO" +
                            " -recalFile %s",
                    1, // just one output file
                    Arrays.asList(md5));
            List<File> result = executeTest("testCountCovariatesNoReadGroups", spec).getFirst();
            paramsFilesNoReadGroupTest.put(bam, result.get(0).getAbsolutePath());
        }
    }

    @Test
    public void testTableRecalibratorNoReadGroups() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12762.SOLID.SRP000031.2009_07.chr1.10_20mb.bam", "32ad300e8c094ed2c1ec6c531180fe70" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFilesNoReadGroupTest.get(bam);
            System.out.printf("PARAMS FOR %s is %s%n", bam, paramsFile);
            if ( paramsFile != null ) {
                WalkerTestSpec spec = new WalkerTestSpec(
                        "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                                " -T TableRecalibration" +
                                " -I " + bam +
                                " -L 1:10,100,000-10,300,000" +
                                " -outputBam %s" +
                                " --no_pg_tag" +
                                " --solid_recal_mode SET_Q_ZERO" +
                                " --default_platform illumina" +
                                " -recalFile " + paramsFile,
                        1, // just one output file
                        Arrays.asList(md5));
                executeTest("testTableRecalibratorNoReadGroups", spec);
            }
        }
    }

    @Test
    public void testCountCovariatesNoIndex() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.noindex.bam", "850f2a2d5bc94cc22b3b038b424252c6" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -cov ReadGroupCovariate" +
                            " -cov QualityScoreCovariate" +
                            " --solid_recal_mode DO_NOTHING" +
                            " -recalFile %s" +
                            " -U",
                    1, // just one output file
                    Arrays.asList(md5));
            List<File> result = executeTest("testCountCovariatesNoIndex", spec).getFirst();
            paramsFiles.put(bam, result.get(0).getAbsolutePath());
        }
    }

    @Test
    public void testTableRecalibratorNoIndex() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.noindex.bam", "aa38b04c6b58badabb6b09d590284a2a" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFiles.get(bam);
            System.out.printf("PARAMS FOR %s is %s%n", bam, paramsFile);
            if ( paramsFile != null ) {
                WalkerTestSpec spec = new WalkerTestSpec(
                        "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                                " -T TableRecalibration" +
                                " -I " + bam +
                                " -outputBam %s" +
                                " --no_pg_tag" +
                                " --solid_recal_mode DO_NOTHING" +
                                " -recalFile " + paramsFile +
                                " -U",
                        1, // just one output file
                        Arrays.asList(md5));
                executeTest("testTableRecalibratorNoIndex", spec);
            }
        }
    }

}
