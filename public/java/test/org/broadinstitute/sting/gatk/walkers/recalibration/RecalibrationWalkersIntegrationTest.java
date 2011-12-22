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
    static HashMap<String, String> paramsFilesSolidIndels = new HashMap<String, String>();

    private static final class CCTest extends TestDataProvider {
        String file, md5;

        private CCTest(final String file, final String md5) {
            super(CCTest.class);
            this.file = file;
            this.md5 = md5;
        }

        public String toString() {
            return "CCTest: " + file;
        }
    }

    @DataProvider(name = "cctestdata")
    public Object[][] createCCTestData() {

        new CCTest( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "9469d6b65880abe4e5babc1c1a69889d" );
        new CCTest( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "17d4b8001c982a70185e344929cf3941");
        new CCTest( validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "36c0c467b6245c2c6c4e9c956443a154" );
        new CCTest( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "ed15f8bf03bb2ea9b7c26844be829c0d" );
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
                        " -knownSites " + b36dbSNP129 +
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

        public String toString() {
            return "TRTest: " + file;
        }
    }

    @DataProvider(name = "trtestdata")
    public Object[][] createTRTestData() {
        new TRTest( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "f020725d9f75ad8f1c14bfae056e250f" );
        new TRTest( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "d04cf1f6df486e45226ebfbf93a188a5");
        new TRTest( validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "b2f4757bc47cf23bd9a09f756c250787" );
        new TRTest( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "313a21a8a88e3460b6e71ec5ffc50f0f" );

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
        e.put( validationDataLocation + "originalQuals.1kg.chr1.1-1K.bam", "bd8288b1fc7629e2e8c2cf7f65fefa8f");

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
                            " -knownSites " + b36dbSNP129,
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testCountCovariatesUseOriginalQuals", spec);
        }
    }

    @Test
    public void testTableRecalibratorMaxQ70() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "f020725d9f75ad8f1c14bfae056e250f" );

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
        e.put( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "1f643bca090478ba68aac88db835a629" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -knownSites " + b36dbSNP129 +
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
        e.put( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "613fb2bbe01af8cbe6a188a10c1582ca" );

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
    public void testCountCovariatesBED() {
        HashMap<String, String> e = new HashMap<String, String>();
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "b00e99219aeafe2516c6232b7d6a0a00");

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -knownSites:bed " + validationDataLocation + "recalibrationTest.bed" +
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
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "7b92788ce92f49415af3a75a2e4a2b33");

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -knownSites:anyNameABCD,VCF3 " + validationDataLocation + "vcfexample3.vcf" +
                            " -T CountCovariates" +
                            " -I " + bam +
                            " -knownSites " + b36dbSNP129 +
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
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.noindex.bam", "f34f7141351a5dbf9664c67260f94e96" );

        for ( Map.Entry<String, String> entry : e.entrySet() ) {
            String bam = entry.getKey();
            String md5 = entry.getValue();

            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + b36KGReference +
                            " -knownSites " + b36dbSNP129 +
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
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.noindex.bam", "13c83656567cee9e93bda9874ee80234" );

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
