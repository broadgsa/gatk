package org.broadinstitute.sting;

import org.apache.commons.io.FileUtils;
import org.apache.log4j.*;
import org.apache.log4j.spi.LoggingEvent;
import org.broadinstitute.sting.commandline.CommandLineUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.testng.Assert;

import javax.swing.*;
import java.io.*;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * User: aaron
 * Date: Apr 14, 2009
 * Time: 10:24:30 AM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date Apr 14, 2009
 * <p/>
 * Class BaseTest
 * <p/>
 * This is the base test class for all of our test cases.  All test cases should extend from this
 * class; it sets up the logger, and resolves the location of directories that we rely on.
 */
@SuppressWarnings("unchecked")
public abstract class BaseTest {
    /** our log, which we want to capture anything from org.broadinstitute.sting */
    public static final Logger logger = CommandLineUtils.getStingLogger();

    public static final String hg18Reference = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta";
    public static final String hg19Reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
    public static final String b36KGReference = "/humgen/1kg/reference/human_b36_both.fasta";
    public static final String b37KGReference = "/humgen/1kg/reference/human_g1k_v37.fasta";
    public static final String GATKDataLocation = "/humgen/gsa-hpprojects/GATK/data/";
    public static final String validationDataLocation = GATKDataLocation + "Validation_Data/";
    public static final String evaluationDataLocation = GATKDataLocation + "Evaluation_Data/";
    public static final String comparisonDataLocation = GATKDataLocation + "Comparisons/";
    public static final String annotationDataLocation = GATKDataLocation + "Annotations/";

    public static final String refseqAnnotationLocation = annotationDataLocation + "refseq/";
    public static final String hg18Refseq = refseqAnnotationLocation + "refGene-big-table-hg18.txt";
    public static final String hg19Refseq = refseqAnnotationLocation + "refGene-big-table-hg19.txt";
    public static final String b36Refseq = refseqAnnotationLocation + "refGene-big-table-b36.txt";
    public static final String b37Refseq = refseqAnnotationLocation + "refGene-big-table-b37.txt";

    public static final String dbsnpDataLocation = GATKDataLocation;
    public static final String hg18dbSNP129 = dbsnpDataLocation + "dbsnp_129_hg18.rod";
    public static final String b36dbSNP129 = dbsnpDataLocation + "dbsnp_129_b36.rod";
    public static final String b37dbSNP129 = dbsnpDataLocation + "dbsnp_129_b37.rod";
    public static final String b37dbSNP132 = dbsnpDataLocation + "dbsnp_132_b37.vcf";

    public static final String hapmapDataLocation = comparisonDataLocation + "Validated/HapMap/3.3/";
    public static final String b37hapmapGenotypes = hapmapDataLocation + "genotypes_r27_nr.b37_fwd.vcf";
    public static final String b37hapmapSites = hapmapDataLocation + "sites_r27_nr.b37_fwd.vcf";

    public static final String intervalsLocation = GATKDataLocation;
    public static final String hg19Intervals = intervalsLocation + "whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list";
    public static final String hg19Chr20Intervals = intervalsLocation + "whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.chr20.interval_list";

    public static final String networkTempDir = "/broad/shptmp/";
    public static final File networkTempDirFile = new File(networkTempDir);

    /**
     * Subdirectory under the ant build directory where we store integration test md5 results
     */
    public static final String MD5_FILE_DB_SUBDIR = "integrationtests";

    public static final String testDir = "public/testdata/";

    /** before the class starts up */
    static {
        // setup a basic log configuration
        CommandLineUtils.configureConsoleLogging();

        // setup our log layout
        PatternLayout layout = new PatternLayout();
        layout.setConversionPattern("TEST %C{1}.%M - %d{HH:mm:ss,SSS} - %m%n");

        // now set the layout of all the loggers to our layout
        CommandLineUtils.setLayout(logger, layout);

        // Set the Root logger to only output warnings.
        logger.setLevel(Level.WARN);

        // find our file sources
        if (!fileExist(hg18Reference) || !fileExist(hg19Reference) || !fileExist(b36KGReference)) {
            logger.fatal("We can't locate the reference directories.  Aborting!");
            throw new RuntimeException("BaseTest setup failed: unable to locate the reference directories");
        }
    }

    /**
     * Simple generic utility class to creating TestNG data providers:
     *
     * 1: inherit this class, as in
     *
     *      private class SummarizeDifferenceTest extends TestDataProvider {
     *         public SummarizeDifferenceTest() {
     *           super(SummarizeDifferenceTest.class);
     *         }
     *         ...
     *      }
     *
     * Provide a reference to your class to the TestDataProvider constructor.
     *
     * 2: Create instances of your subclass.  Return from it the call to getTests, providing
     * the class type of your test
     *
     * @DataProvider(name = "summaries")
     * public Object[][] createSummaries() {
     *   new SummarizeDifferenceTest().addDiff("A", "A").addSummary("A:2");
     *   new SummarizeDifferenceTest().addDiff("A", "B").addSummary("A:1", "B:1");
     *   return SummarizeDifferenceTest.getTests(SummarizeDifferenceTest.class);
     * }
     *
     * This class magically tracks created objects of this
     */
    public static class TestDataProvider {
        private static final Map<Class, List<Object>> tests = new HashMap<Class, List<Object>>();

        /**
         * Create a new TestDataProvider instance bound to the class variable C
         * @param c
         */
        public TestDataProvider(Class c) {
            if ( ! tests.containsKey(c) )
                tests.put(c, new ArrayList<Object>());
            tests.get(c).add(this);
        }

        /**
         * Return all of the data providers in the form expected by TestNG of type class C
         * @param c
         * @return
         */
        public static Object[][] getTests(Class c) {
            List<Object[]> params2 = new ArrayList<Object[]>();
            for ( Object x : tests.get(c) ) params2.add(new Object[]{x});
            return params2.toArray(new Object[][]{});
        }
    }

    /**
     * test if the file exists
     *
     * @param file name as a string
     * @return true if it exists
     */
    public static boolean fileExist(String file) {
        File temp = new File(file);
        return temp.exists();
    }
    
    /**
     * this appender looks for a specific message in the log4j stream.
     * It can be used to verify that a specific message was generated to the logging system.
     */
    public static class ValidationAppender extends AppenderSkeleton {

        private boolean foundString = false;
        private String targetString = "";

        public ValidationAppender(String target) {
            targetString = target;
        }

        @Override
        protected void append(LoggingEvent loggingEvent) {
            if (loggingEvent.getMessage().equals(targetString))
                foundString = true;
        }

        public void close() {
            // do nothing
        }

        public boolean requiresLayout() {
            return false;
        }

        public boolean foundString() {
            return foundString;
        }
    }

    /**
     * a little utility function for all tests to md5sum a file
     * Shameless taken from:
     *
     * http://www.javalobby.org/java/forums/t84420.html
     *
     * @param file the file
     * @return a string
     */
    public static String md5SumFile(File file) {
        MessageDigest digest;
        try {
            digest = MessageDigest.getInstance("MD5");
        } catch (NoSuchAlgorithmException e) {
            throw new ReviewedStingException("Unable to find MD5 digest");
        }
        InputStream is;
        try {
            is = new FileInputStream(file);
        } catch (FileNotFoundException e) {
            throw new ReviewedStingException("Unable to open file " + file);
        }
        byte[] buffer = new byte[8192];
        int read;
        try {
            while ((read = is.read(buffer)) > 0) {
                digest.update(buffer, 0, read);
            }
            byte[] md5sum = digest.digest();
            BigInteger bigInt = new BigInteger(1, md5sum);
            return bigInt.toString(16);

        }
        catch (IOException e) {
            throw new ReviewedStingException("Unable to process file for MD5", e);
        }
        finally {
            try {
                is.close();
            }
            catch (IOException e) {
                throw new ReviewedStingException("Unable to close input stream for MD5 calculation", e);
            }
        }
    }

    protected static void ensureMd5DbDirectory() {
        // todo -- make path
        File dir = new File(MD5_FILE_DB_SUBDIR);
        if ( ! dir.exists() ) {
            System.out.printf("##### Creating MD5 db %s%n", MD5_FILE_DB_SUBDIR);
            if ( ! dir.mkdir() ) {
                throw new ReviewedStingException("Infrastructure failure: failed to create md5 directory " + MD5_FILE_DB_SUBDIR);
            }
        }
    }

    protected static File getFileForMD5(final String md5) {
        final String basename = String.format("%s.integrationtest", md5);
        return new File(MD5_FILE_DB_SUBDIR + "/" + basename);
    }

    private static void updateMD5Db(final String md5, final File resultsFile) {
        // todo -- copy results file to DB dir if needed under filename for md5
        final File dbFile = getFileForMD5(md5);
        if ( ! dbFile.exists() ) {
            // the file isn't already in the db, copy it over
            System.out.printf("##### Updating MD5 file: %s%n", dbFile.getPath());
            try {
                FileUtils.copyFile(resultsFile, dbFile);
            } catch ( IOException e ) {
                throw new ReviewedStingException(e.getMessage());
            }
        } else {
            System.out.printf("##### MD5 file is up to date: %s%n", dbFile.getPath());

        }
    }

    private static String getMD5Path(final String md5, final String valueIfNotFound) {
        // todo -- look up the result in the directory and return the path if it exists
        final File dbFile = getFileForMD5(md5);
        return dbFile.exists() ? dbFile.getPath() : valueIfNotFound;
    }

    public static byte[] getBytesFromFile(File file) throws IOException {
        InputStream is = new FileInputStream(file);

        // Get the size of the file
        long length = file.length();

        if (length > Integer.MAX_VALUE) {
            // File is too large
        }

        // Create the byte array to hold the data
        byte[] bytes = new byte[(int) length];

        // Read in the bytes
        int offset = 0;
        int numRead = 0;
        while (offset < bytes.length
                && (numRead = is.read(bytes, offset, bytes.length - offset)) >= 0) {
            offset += numRead;
        }

        // Ensure all the bytes have been read in
        if (offset < bytes.length) {
            throw new IOException("Could not completely read file " + file.getName());
        }

        // Close the input stream and return bytes
        is.close();
        return bytes;
    }

    /**
     * Tests a file MD5 against an expected value, returning the MD5.  NOTE: This function WILL throw an exception if the MD5s are different.
     * @param name Name of the test.
     * @param resultsFile File to MD5.
     * @param expectedMD5 Expected MD5 value.
     * @param parameterize If true or if expectedMD5 is an empty string, will print out the calculated MD5 instead of error text.
     * @return The calculated MD5.
     */
    public static String assertMatchingMD5(final String name, final File resultsFile, final String expectedMD5, final boolean parameterize) {
        String filemd5sum = testFileMD5(name, resultsFile, expectedMD5, parameterize);
        
        if (parameterize || expectedMD5.equals("")) {
            // Don't assert
        } else {
            Assert.assertEquals(filemd5sum, expectedMD5, name + " Mismatching MD5s");
            System.out.println(String.format("  => %s PASSED", name));
        }

        return filemd5sum;
    }


    /**
     * Tests a file MD5 against an expected value, returning the MD5.  NOTE: This function WILL NOT throw an exception if the MD5s are different.
     * @param name Name of the test.
     * @param resultsFile File to MD5.
     * @param expectedMD5 Expected MD5 value.
     * @param parameterize If true or if expectedMD5 is an empty string, will print out the calculated MD5 instead of error text.
     * @return The calculated MD5.
     */
    public static String testFileMD5(final String name, final File resultsFile, final String expectedMD5, final boolean parameterize) {
        try {
            byte[] bytesOfMessage = getBytesFromFile(resultsFile);
            byte[] thedigest = MessageDigest.getInstance("MD5").digest(bytesOfMessage);
            BigInteger bigInt = new BigInteger(1, thedigest);
            String filemd5sum = bigInt.toString(16);
            while (filemd5sum.length() < 32) filemd5sum = "0" + filemd5sum; // pad to length 32

            //
            // copy md5 to integrationtests
            //
            updateMD5Db(filemd5sum, resultsFile);

            if (parameterize || expectedMD5.equals("")) {
                System.out.println(String.format("PARAMETERIZATION[%s]: file %s has md5 = %s, stated expectation is %s, equal? = %b",
                                                 name, resultsFile, filemd5sum, expectedMD5, filemd5sum.equals(expectedMD5)));
            } else {
                System.out.println(String.format("Checking MD5 for %s [calculated=%s, expected=%s]", resultsFile, filemd5sum, expectedMD5));
                System.out.flush();

                if ( ! expectedMD5.equals(filemd5sum) ) {
                    // we are going to fail for real in assertEquals (so we are counted by the testing framework).
                    // prepare ourselves for the comparison
                    System.out.printf("##### Test %s is going fail #####%n", name);
                    String pathToExpectedMD5File = getMD5Path(expectedMD5, "[No DB file found]");
                    String pathToFileMD5File = getMD5Path(filemd5sum, "[No DB file found]");
                    System.out.printf("##### Path to expected   file (MD5=%s): %s%n", expectedMD5, pathToExpectedMD5File);
                    System.out.printf("##### Path to calculated file (MD5=%s): %s%n", filemd5sum, pathToFileMD5File);
                    System.out.printf("##### Diff command: diff %s %s%n", pathToExpectedMD5File, pathToFileMD5File);

                    // todo -- add support for simple inline display of the first N differences for text file
                }
            }

            return filemd5sum;
        } catch (Exception e) {
            throw new RuntimeException("Failed to read bytes from calls file: " + resultsFile, e);
        }
    }

    /**
     * Creates a temp file that will be deleted on exit after tests are complete.
     * @param name Prefix of the file.
     * @param extension Extension to concat to the end of the file.
     * @return A file in the temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     */
    public static File createTempFile(String name, String extension) {
        try {
            File file = File.createTempFile(name, extension);
            file.deleteOnExit();
            return file;
        } catch (IOException ex) {
            throw new ReviewedStingException("Cannot create temp file: " + ex.getMessage(), ex);
        }
    }

    /**
     * Creates a temp file that will be deleted on exit after tests are complete.
     * @param name Prefix of the file.
     * @param extension Extension to concat to the end of the file.
     * @return A file in the network temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     */
    public static File createNetworkTempFile(String name, String extension) {
        try {
            File file = File.createTempFile(name, extension, networkTempDirFile);
            file.deleteOnExit();
            return file;
        } catch (IOException ex) {
            throw new ReviewedStingException("Cannot create temp file: " + ex.getMessage(), ex);
        }
    }
}
