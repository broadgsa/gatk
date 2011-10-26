package org.broadinstitute.sting;

import org.apache.commons.io.FileUtils;
import org.apache.log4j.*;
import org.apache.log4j.spi.LoggingEvent;
import org.broadinstitute.sting.commandline.CommandLineUtils;
import org.broadinstitute.sting.gatk.walkers.diffengine.DiffEngine;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.testng.Assert;

import javax.swing.*;
import java.io.*;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;

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
    //public static final String b37KGReference = "/Users/depristo/Desktop/broadLocal/localData/human_g1k_v37.fasta";
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
    public static final String b36dbSNP129 = dbsnpDataLocation + "dbsnp_129_b36.vcf";
    public static final String b37dbSNP129 = dbsnpDataLocation + "dbsnp_129_b37.vcf";
    public static final String b37dbSNP132 = dbsnpDataLocation + "dbsnp_132_b37.vcf";
    public static final String hg18dbSNP132 = dbsnpDataLocation + "dbsnp_132.hg18.vcf";

    public static final String hapmapDataLocation = comparisonDataLocation + "Validated/HapMap/3.3/";
    public static final String b37hapmapGenotypes = hapmapDataLocation + "genotypes_r27_nr.b37_fwd.vcf";
    public static final String b37hapmapSites = hapmapDataLocation + "sites_r27_nr.b37_fwd.vcf";

    public static final String intervalsLocation = GATKDataLocation;
    public static final String hg19Intervals = intervalsLocation + "whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list";
    public static final String hg19Chr20Intervals = intervalsLocation + "whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.chr20.interval_list";

    public static final String networkTempDir = "/broad/shptmp/";
    public static final File networkTempDirFile = new File(networkTempDir);

    public static final File testDirFile = new File("public/testdata/");
    public static final String testDir = testDirFile.getAbsolutePath() + "/";

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
//        if (!fileExist(hg18Reference) || !fileExist(hg19Reference) || !fileExist(b36KGReference)) {
//            logger.fatal("We can't locate the reference directories.  Aborting!");
//            throw new RuntimeException("BaseTest setup failed: unable to locate the reference directories");
//        }
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
        private final String name;

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
