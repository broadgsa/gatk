package org.broadinstitute.sting;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import org.apache.log4j.AppenderSkeleton;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.apache.log4j.spi.LoggingEvent;
import org.broadinstitute.sting.commandline.CommandLineUtils;
import org.broadinstitute.sting.utils.crypt.CryptUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.io.IOUtils;

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

            
    public static String hg18Reference;
    public static String hg19Reference;
    public static String b36KGReference;
    public static String b37KGReference;
    public static String GATKDataLocation;
    public static String validationDataLocation;
    public static String evaluationDataLocation;
    public static String comparisonDataLocation;
    public static String annotationDataLocation;

    public static String b37GoodBAM;
    public static String b37GoodNA12878BAM;
    public static String b37_NA12878_OMNI;

    public static String refseqAnnotationLocation;
    public static String hg18Refseq;
    public static String hg19Refseq;
    public static String b36Refseq;
    public static String b37Refseq;

    public static String dbsnpDataLocation;
    public static String b36dbSNP129;
    public static String b37dbSNP129;
    public static String b37dbSNP132;
    public static String hg18dbSNP132;

    public static String hapmapDataLocation;
    public static String b37hapmapGenotypes;
    public static String b37hapmapSites;

    public static String intervalsLocation;
    public static String hg19Intervals;
    public static String hg19Chr20Intervals;

    public static boolean REQUIRE_NETWORK_CONNECTION;
    public static String networkTempDirRoot;
    public static String networkTempDir;
    public static File networkTempDirFile;

    public static final File testDirFile = new File("public/testdata/");
    public static final String testDir = testDirFile.getAbsolutePath() + "/";

    public static final String keysDataLocation = validationDataLocation + "keys/";
    public static final String gatkKeyFile = CryptUtils.GATK_USER_KEY_DIRECTORY + "gsamembers_broadinstitute.org.key";

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
        
        // Get and set all file paths from the config file.
        setupResourcePaths();
        
        if ( REQUIRE_NETWORK_CONNECTION ) {        
            networkTempDirFile = IOUtils.tempDir("temp.", ".dir", new File(networkTempDirRoot + System.getProperty("user.name")));
            networkTempDirFile.deleteOnExit();
            networkTempDir = networkTempDirFile.getAbsolutePath() + "/";

            // find our file sources
            if (!fileExist(hg18Reference) || !fileExist(hg19Reference) || !fileExist(b36KGReference)) {
                logger.fatal("We can't locate the reference directories.  Aborting!");
                throw new RuntimeException("BaseTest setup failed: unable to locate the reference directories");
            }
        } else {
            networkTempDir = null;
            networkTempDirFile = null;
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
     * @DataProvider(name = "summaries"
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
        protected String name;

        /**
         * Create a new TestDataProvider instance bound to the class variable C
         * @param c
         */
        public TestDataProvider(Class c, String name) {
            if ( ! tests.containsKey(c) )
                tests.put(c, new ArrayList<Object>());
            tests.get(c).add(this);
            this.name = name;
        }

        public TestDataProvider(Class c) {
            this(c, "");
        }

        public void setName(final String name) {
            this.name = name;
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

        @Override
        public String toString() {
            return "TestDataProvider("+name+")";
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
     * Setting up paths to all resources from the external test.conf file.
     */
    private static void setupResourcePaths() {
    /** Load settings from external configuration file */
    Properties properties = new Properties();

    FileInputStream configFile;
    try {
        /** Get the settings file from the test.conf file in settings under the gatk root. */
        configFile = new FileInputStream(new File("settings/test.conf"));
        properties.load(configFile);
        configFile.close();

        hg18Reference = properties.getProperty("hg18Reference");
        hg19Reference = properties.getProperty("hg19Reference");
        b36KGReference = properties.getProperty("b36KGReference");	       
        b37KGReference = properties.getProperty("b37KGReference");
        GATKDataLocation = properties.getProperty("GATKDataLocation");
        validationDataLocation = properties.getProperty("validationDataLocation");
        evaluationDataLocation = properties.getProperty("evaluationDataLocation");
        comparisonDataLocation = properties.getProperty("comparisonDataLocation");
        annotationDataLocation = properties.getProperty("annotationDataLocation");

        b37GoodBAM = properties.getProperty("b37GoodBAM");
        b37GoodNA12878BAM = properties.getProperty("b37GoodNA12878BAM");
        b37_NA12878_OMNI = properties.getProperty("b37_NA12878_OMNI");

        refseqAnnotationLocation = properties.getProperty("refseqAnnotationLocation");
        hg18Refseq = properties.getProperty("hg18Refseq");
        hg19Refseq = properties.getProperty("hg19Refseq");
        b36Refseq = properties.getProperty("b36Refseq");
        b37Refseq = properties.getProperty("b37Refseq");

        dbsnpDataLocation = properties.getProperty("dbsnpDataLocation");
        b36dbSNP129 = properties.getProperty("b36dbSNP129");
        b37dbSNP129 = properties.getProperty("b37dbSNP129");
        b37dbSNP132 = properties.getProperty("b37dbSNP132");
        hg18dbSNP132 = properties.getProperty("hg18dbSNP132");

        hapmapDataLocation = properties.getProperty("hapmapDataLocatio");
        b37hapmapGenotypes = properties.getProperty("b37hapmapGenotypes");
        b37hapmapSites = properties.getProperty("b37hapmapSites");

        intervalsLocation = properties.getProperty("intervalsLocation");
        hg19Intervals = properties.getProperty("hg19Intervals");
        hg19Chr20Intervals = properties.getProperty("hg19Chr20Intervals");

        REQUIRE_NETWORK_CONNECTION = Boolean.getBoolean(properties.getProperty("REQUIRE_NETWORK_CONNECTION"));
        networkTempDirRoot = properties.getProperty("networkTempDirRoot");

    }
    catch (FileNotFoundException e) {
        logger.fatal("Could not find settings file \"test.conf\" in the settings directory.");
        throw new RuntimeException("BaseTest setup failed: could not find \"test.conf\" file in the settings directory.");
        } 
    catch (IOException e) {
        logger.fatal("Could not read \"test.conf\" file. Aborting!");
        throw new RuntimeException("BaseTest setup failed: tried to get settings from test.conf, but could not read it.");
        }		
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
     * @param name Name of the file.
     * @return A file in the network temporary directory with name, which will be deleted after the program exits.
     */
    public static File createNetworkTempFile(String name) {
        File file = new File(networkTempDirFile, name);
        file.deleteOnExit();
        return file;
    }
}
