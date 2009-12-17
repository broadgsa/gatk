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
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12892.SLX.SRP000031.2009_06.selected.bam", "16f87c9644b27c3c3fe7a963eed45d2d" );
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "5a4873064e602e1f945f3881c6b4a3e5");
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "b039f405580ad844cd1dfd4d8ba25913" );
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "dbbb9c49f2e0b1f117b88eb386fb6d0f" );

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
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12892.SLX.SRP000031.2009_06.selected.bam", "d7f3e0db5ed9fefc917144a0f937d50d" );
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "1f1770363f381d06e8db9fead03ad7bb");
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "3fe6ee95cfa6155f8cc8e3e0a6330ad2" );
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "732267e9cd02a72b8c85238f47487509" );

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
                                " --solid_recal_mode SET_Q_ZERO" +
                                " -recalFile " + paramsFile,
                        1, // just one output file
                        Arrays.asList(md5));
                executeTest("testTableRecalibrator1", spec);
            }
        }
    }

    @Test
    public void testCountCovariatesVCF() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "bc0e62dc366f2f6ca70847a832b4aad4");

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R /broad/1KG/reference/human_b36_both.fasta" +
                            " -B dbsnp,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/vcfexample3.vcf" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -L 1:10,000,000-10,200,000" +
                            " -cov ReadGroupCovariate" +
                            " -cov QualityScoreCovariate" +
                            " -cov CycleCovariate" +
                            " -cov DinucCovariate" +
                            " --sorted_output" +
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
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12762.SOLID.SRP000031.2009_07.chr1.10_20mb.bam", "9d83a1653e076e028fecf4c2575f56e9" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R /broad/1KG/reference/human_b36_both.fasta" +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -L 1:10,000,000-10,200,000" +
                            " -cov ReadGroupCovariate" +
                            " -cov QualityScoreCovariate" +
                            " -cov CycleCovariate" +
                            " -cov DinucCovariate" +
                            " --default_platform illumina" +
                            " --sorted_output" +
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
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12762.SOLID.SRP000031.2009_07.chr1.10_20mb.bam", "8853add9cfd0a8171ca917b1c6452e0a" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFilesNoReadGroupTest.get(bam);
            System.out.printf("PARAMS FOR %s is %s%n", bam, paramsFile);
            if ( paramsFile != null ) {
                WalkerTestSpec spec = new WalkerTestSpec(
                        "-R /broad/1KG/reference/human_b36_both.fasta" +
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
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.noindex.bam", "f5a0adce5458ab661029fa5056a55551" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R /broad/1KG/reference/human_b36_both.fasta" +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -cov ReadGroupCovariate" +
                            " -cov QualityScoreCovariate" +
                            " --sorted_output" +
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
        e.put( "/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.noindex.bam", "dd8705b54f0ffbedb689608221ade13a" );

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
