package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
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

    private static final class CCTest extends TestDataProvider {
        String file, md5;

        private CCTest(final String file, final String md5) {
            super(CCTest.class);
            this.file = file;
            this.md5 = md5;
        }
    }

    @DataProvider(name = "cctestdata")
    public Object[][] createCCTestData() {

        new CCTest( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "" );
        new CCTest( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "");
        new CCTest( validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "" );
        new CCTest( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "" );
        return CCTest.getTests(CCTest.class);
    }

    @Test(dataProvider = "cctestdata")
    public void testCountCovariates1(CCTest test) {
        testCC(test, "");
    }

    @Test(dataProvider = "cctestdata")
    public void testCountCovariates4(CCTest test) {
        testCC(test, " -nt 4");
    }

    private final void testCC(CCTest test, String parallelism) {
        String bam = test.file;
        String md5 = test.md5;

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -B:dbsnp,vcf " + GATKDataLocation + "dbsnp_132.b36.excluding_sites_after_129.vcf" +
                        " -T CountCovariates" +
                        " -I " + bam +
                        ( bam.equals( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" )
                                ? " -L 1:10,800,000-10,810,000" : " -L 1:10,000,000-10,200,000" ) +
                        " -cov ReadGroupCovariate" +
                        " -cov QualityScoreCovariate" +
                        " -cov CycleCovariate" +
                        " -cov DinucCovariate" +
                        " --solid_recal_mode SET_Q_ZERO" +
                        " -recalFile %s" + parallelism,
                1, // just one output file
                Arrays.asList(md5));
        List<File> result = executeTest("testCountCovariates1" + parallelism, spec).getFirst();
        paramsFiles.put(bam, result.get(0).getAbsolutePath());
    }


    private static final class TRTest extends TestDataProvider {
        String file, md5;

        private TRTest(final String file, final String md5) {
            super(TRTest.class);
            this.file = file;
            this.md5 = md5;
        }
    }

    @DataProvider(name = "trtestdata")
    public Object[][] createTRTestData() {
        new TRTest( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "" );
        new TRTest( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "");
        new TRTest( validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "" );
        new TRTest( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "" );
        return TRTest.getTests(TRTest.class);
    }

    @Test(dataProvider = "trtestdata", dependsOnMethods = "testCountCovariates1")
    public void testTableRecalibrator1(TRTest test) {
        String bam = test.file;
        String md5 = test.md5;
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

    @Test
    public void testCountCovariatesUseOriginalQuals() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "originalQuals.1kg.chr1.1-1K.bam", "");

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
                            " -B:dbsnp,vcf " + GATKDataLocation + "dbsnp_132.b36.excluding_sites_after_129.vcf",
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testCountCovariatesUseOriginalQuals", spec);
        }
    }

    @Test
    public void testTableRecalibratorMaxQ70() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "" );

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
        e.put( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -B:dbsnp,vcf " + GATKDataLocation + "dbsnp_132.b36.excluding_sites_after_129.vcf" +
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
        e.put( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "" );

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
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "");

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -B:dbsnp,VCF3 " + validationDataLocation + "vcfexample3.vcf" +
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
    public void testCountCovariatesBED() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "");

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -B:bed,bed " + validationDataLocation + "recalibrationTest.bed" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -L 1:10,000,000-10,200,000" +
                            " -standard" +
                            " --solid_recal_mode SET_Q_ZERO" +
                            " -recalFile %s",
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testCountCovariatesBED", spec);
        }
    }

    @Test
    public void testCountCovariatesVCFPlusDBsnp() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "");

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -B:anyNameABCD,VCF3 " + validationDataLocation + "vcfexample3.vcf" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -B:dbsnp,vcf " + GATKDataLocation + "dbsnp_132.b36.excluding_sites_after_129.vcf" +
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
    public void testCountCovariatesNoIndex() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.noindex.bam", "" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -B:dbsnp,vcf " + GATKDataLocation + "dbsnp_132.b36.excluding_sites_after_129.vcf" +
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
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.noindex.bam", "" );

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
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "");

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
