package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;
import java.util.List;
import java.io.File;

public class RecalibrationWalkersIntegrationTest extends WalkerTest {
    static HashMap<String, String> paramsFiles = new HashMap<String, String>();
    static HashMap<String, String> paramsFilesNoReadGroupTest = new HashMap<String, String>();
    static HashMap<String, String> paramsFilesSolidIndels = new HashMap<String, String>();

    @Test
    public void testCountCovariates1() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "e2ea7507feb66651f52320c5a46433c2" );
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "b2191ea11f528b9605b727d8a73dd1e1");
        e.put( validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "596a9ec9cbc1da70481e45a5a588a41d" );
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "507dbd3ba6f54e066d04c4d24f59c3ab" );

        for ( String parallelism : Arrays.asList("", " -nt 4")) {
            for ( Map.Entry<String, String> entry : e.entrySet() ) {
                String bam = entry.getKey();
                String md5 = entry.getValue();

                WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                        "-R " + b36KGReference +
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
                                " -recalFile %s" + parallelism,
                        1, // just one output file
                        Arrays.asList(md5));
                List<File> result = executeTest("testCountCovariates1" + parallelism, spec).getFirst();
                paramsFiles.put(bam, result.get(0).getAbsolutePath());
            }
        }
    }
    
    @Test
    public void testTableRecalibrator1() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "da79dcc74336a6333f75aeb231590436" );
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "6c5ac8f71a77f09eec3e38932ba274cf");
        e.put( validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "b26a07cc3d3f6a607e786ebcc2c9c8ac" );
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "504152c1a94d6998af9785bcbddbb3d6" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFiles.get(bam);
            System.out.printf("PARAMS FOR %s is %s%n", bam, paramsFile);
            if ( paramsFile != null ) {
                WalkerTestSpec spec = new WalkerTestSpec(
                        "-R " + b36KGReference +
                                " -T TableRecalibration" +
                                " -I " + bam +
                                ( bam.equals( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" )
                                    ? " -L 1:10,800,000-10,810,000" : " -L 1:10,100,000-10,300,000" ) +
                                " -o %s" +
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
    public void testCountCovariatesUseOriginalQuals() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "originalQuals.1kg.chr1.1-1K.bam", "db6c0dca1ec121f8c2e6d20f4969dee5");

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -L 1:1-1,000" +
                            " -standard" +
                            " -OQ" +
                            " -recalFile %s" +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod",
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testCountCovariatesUseOriginalQuals", spec);
        }
    }

    @Test
    public void testTableRecalibratorMaxQ70() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "4d24b4683a9b20e93837da75fff2489f" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFiles.get(bam);
            System.out.printf("PARAMS FOR %s is %s%n", bam, paramsFile);
            if ( paramsFile != null ) {
                WalkerTestSpec spec = new WalkerTestSpec(
                        "-R " + b36KGReference +
                                " -T TableRecalibration" +
                                " -I " + bam +
                                ( bam.equals( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" )
                                    ? " -L 1:10,800,000-10,810,000" : " -L 1:10,100,000-10,300,000" ) +
                                " -o %s" +
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
    public void testCountCovariatesSolidIndelsRemoveRefBias() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "013822cfa4f276d48ca99c014c23c124" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -standard" +
                            " -U" +
                            " -L 1:10,000,000-20,000,000" +
                            " --solid_recal_mode REMOVE_REF_BIAS" +
                            " -recalFile %s",
                    1, // just one output file
                    Arrays.asList(md5));
            List<File> result = executeTest("testCountCovariatesSolidIndelsRemoveRefBias", spec).getFirst();
            paramsFilesSolidIndels.put(bam, result.get(0).getAbsolutePath());
        }
    }

    @Test
    public void testTableRecalibratorSolidIndelsRemoveRefBias() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "85c0500f9f4ba0d613a77d34b5c84a64" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFilesSolidIndels.get(bam);
            System.out.printf("PARAMS FOR %s is %s%n", bam, paramsFile);
            if ( paramsFile != null ) {
                WalkerTestSpec spec = new WalkerTestSpec(
                        "-R " + b36KGReference +
                                " -T TableRecalibration" +
                                " -I " + bam +
                                " -o %s" +
                                " --no_pg_tag" +
                                " -U" +
                                " -L 1:10,000,000-20,000,000" +
                                " --solid_recal_mode REMOVE_REF_BIAS" +
                                " -recalFile " + paramsFile,
                        1, // just one output file
                        Arrays.asList(md5));
                executeTest("testTableRecalibratorSolidIndelsRemoveRefBias", spec);
            }
        }
    }




    @Test
    public void testCountCovariatesVCF() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "3700eaf567e4937f442fc777a226d6ad");

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -B:dbsnp,VCF " + validationDataLocation + "vcfexample3.vcf" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -L 1:10,000,000-10,200,000" +
                            " -standard" +
                            " --solid_recal_mode SET_Q_ZERO" +
                            " -recalFile %s",
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testCountCovariatesVCF", spec);
        }
    }

    @Test
    public void testCountCovariatesVCFPlusDBsnp() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "949c2ecb24a4189e106d372b05ec725f");

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -B:anyNameABCD,VCF " + validationDataLocation + "vcfexample3.vcf" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -L 1:10,000,000-10,200,000" +
                            " -cov ReadGroupCovariate" +
                            " -cov QualityScoreCovariate" +
                            " -cov CycleCovariate" +
                            " -cov DinucCovariate" +
                            " --solid_recal_mode SET_Q_ZERO" +
                            " -recalFile %s",
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testCountCovariatesVCFPlusDBsnp", spec);
        }
    }

    @Test
    public void testCountCovariatesNoReadGroups() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12762.SOLID.SRP000031.2009_07.chr1.10_20mb.bam", "62dab3db2172695cf95fba7f543a4058" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -L 1:10,000,000-10,200,000" +
                            " -cov ReadGroupCovariate" +
                            " -cov QualityScoreCovariate" +
                            " -cov CycleCovariate" +
                            " -cov DinucCovariate" +
                            " --default_read_group DefaultReadGroup" +
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
        e.put( validationDataLocation + "NA12762.SOLID.SRP000031.2009_07.chr1.10_20mb.bam", "a2d7b780052fe31b0d37fc6462ab58ca" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFilesNoReadGroupTest.get(bam);
            System.out.printf("PARAMS FOR %s is %s%n", bam, paramsFile);
            if ( paramsFile != null ) {
                WalkerTestSpec spec = new WalkerTestSpec(
                        "-R " + b36KGReference +
                                " -T TableRecalibration" +
                                " -I " + bam +
                                " -L 1:10,100,000-10,300,000" +
                                " -o %s" +
                                " --no_pg_tag" +
                                " --solid_recal_mode SET_Q_ZERO" +
                                " --default_read_group DefaultReadGroup" +
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
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.noindex.bam", "abc4248cb5f718594a63409a151d679e" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
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
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.noindex.bam", "4ae0b3d4e2f6a7c1c60ae64944afff5e" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();
            String paramsFile = paramsFiles.get(bam);
            System.out.printf("PARAMS FOR %s is %s%n", bam, paramsFile);
            if ( paramsFile != null ) {
                WalkerTestSpec spec = new WalkerTestSpec(
                        "-R " + b36KGReference +
                                " -T TableRecalibration" +
                                " -I " + bam +
                                " -o %s" +
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


    @Test
    public void testCountCovariatesFailWithoutDBSNP() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SOLID.bam", "");

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -L 1:10,000,000-10,200,000" +
                            " -standard" +
                            " --solid_recal_mode SET_Q_ZERO" +
                            " -recalFile %s",
                    1, // just one output file
                    UserException.CommandLineException.class);
            executeTest("testCountCovariatesFailWithoutDBSNP", spec);
        }
    }

}
