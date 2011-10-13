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
        new CCTest( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "5a52b00d9794d27af723bcf93366681e" );
        new CCTest( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "17d4b8001c982a70185e344929cf3941");
        new CCTest( validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "714e65d6cb51ae32221a77ce84cbbcdc" );
        new CCTest( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "64e9f17a1cf6fc04c1f2717c2d2eca67" );
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
        new TRTest( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "2864f231fab7030377f3c8826796e48f" );
        new TRTest( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "d04cf1f6df486e45226ebfbf93a188a5");
        new TRTest( validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam", "74314e5562c1a65547bb0edaacffe602" );
        new TRTest( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam", "2a37c6001826bfabf87063b1dfcf594f" );
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
        e.put( validationDataLocation + "originalQuals.1kg.chr1.1-1K.bam", "278846c55d97bd9812b758468a83f559");

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
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "2864f231fab7030377f3c8826796e48f" );

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
        e.put( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "8379f24cf5312587a1f92c162ecc220f" );

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
        e.put( validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "2ad4c17ac3ed380071137e4e53a398a5" );

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
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "b460478d9683e827784e42bc352db8bb");

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
        e.put( validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam", "9131d96f39badbf9753653f55b148012");

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
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.noindex.bam", "8993d32df5cb66c7149f59eccbd57f4c" );

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
        e.put( validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.noindex.bam", "5f913c98ca99754902e9d34f99df468f" );

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
