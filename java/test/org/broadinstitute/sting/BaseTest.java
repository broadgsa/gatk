package org.broadinstitute.sting;

import org.apache.log4j.*;
import org.apache.log4j.spi.LoggingEvent;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.junit.*;

import java.io.*;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Enumeration;

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
public abstract class BaseTest {
    /** our log, which we want to capture anything from org.broadinstitute.sting */
    public static Logger logger = Logger.getRootLogger();

    protected static String hg18Reference = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta";
    protected static String hg19Reference = "/seq/references/Homo_sapiens_assembly19/v0/Homo_sapiens_assembly19.fasta";
    protected static String b36KGReference = "/humgen/1kg/reference/human_b36_both.fasta";
    protected static String b37KGReference = "/humgen/1kg/reference/human_g1k_v37.fasta";
    protected static String GATKDataLocation = "/humgen/gsa-hpprojects/GATK/data/";
    protected static String validationDataLocation = GATKDataLocation + "Validation_Data/";
    protected static String evaluationDataLocation = GATKDataLocation + "Evaluation_Data/";
    protected static String comparisonDataLocation = GATKDataLocation + "Comparisons/";

    protected static String testDir = "testdata/";
    protected static boolean alreadySetup = false;
    

    /** before the class starts up */
    @BeforeClass
    public static void baseStartup() {
        if (!alreadySetup) {

            alreadySetup = true;
            // setup a basic log configuration
            BasicConfigurator.configure();

            // setup our log layout
            PatternLayout layout = new PatternLayout();
            layout.setConversionPattern("TEST %C{1}.%M - %d{HH:mm:ss,SSS} - %m%n");

            // now set the layout of all the loggers to our layout
            Enumeration<Appender> en = logger.getAllAppenders();
            for (; en.hasMoreElements();) {
                Appender app = en.nextElement();
                app.setLayout(layout);
            }
            logger.setLevel(Level.WARN);

            // find our file sources
            if (!fileExist(hg18Reference) || !fileExist(hg19Reference) || !fileExist(b36KGReference)) {
                logger.fatal("We can't locate the reference directories.  Aborting!");
                throw new RuntimeException("BaseTest setup failed: unable to locate the reference directories");
            }
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


}
