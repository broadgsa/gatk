package org.broadinstitute.sting;

import org.apache.log4j.*;
import org.apache.log4j.spi.LoggingEvent;
import org.broadinstitute.sting.commandline.CommandLineUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.*;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

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

    public final String testDir = "testdata/";

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
