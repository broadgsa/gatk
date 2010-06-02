package org.broadinstitute.sting;

import org.apache.log4j.*;
import org.apache.log4j.spi.LoggingEvent;
import org.junit.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
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

    protected static String seqLocation = "/seq/";
    protected static String oneKGLocation = "/broad/1KG/";
    protected static String hg18Reference = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta";
    protected static String b36KGReference = "/broad/1KG/reference/human_b36_both.fasta";
    protected static String GATKDataLocation = "/humgen/gsa-hpprojects/GATK/data/";
    protected static String validationDataLocation = GATKDataLocation + "Validation_Data/";
    protected static String evaluationDataLocation = GATKDataLocation + "Evaluation_Data/";

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
            try {
                findFileLocations();
            } catch (IOException e) {
                logger.fatal("We can't locate the base /seq and /broad/1KG directories, for the following reason: " + e.getMessage());
                throw new RuntimeException("BaseTest setup failed: findFileLocations emited exception",e);
            }
        }
    }
    
    public static void findFileLocations() throws IOException {
        // if either doesn't exist
        if (!fileExist(seqLocation) || !fileExist(oneKGLocation)) {

            // get the working directory
            String workDir = System.getProperty("user.dir");

            if (!fileExist("test.conf")) {
                throw new IOException("Unable to find both data directories or your test.conf, make sure it's placed in the base Sting directory");
            }

            BufferedReader inputStream = new BufferedReader(new FileReader("test.conf"));

            String line = "";
            while ((line = inputStream.readLine()) != null) {
                String[] array = line.split("=");

                // check the length
                if (array.length != 2) {
                    throw new IOException("Line : " + line + ", split != 2");
                }
                // clean up the right side
                if (array[1].contains("\n")) {
                    array[1] = array[1].substring(0, array[1].indexOf("\n") - 1);
                }
                // check to see if the right side exists
                if (!fileExist(array[1])) {
                    throw new IOException("Line : " + line + ", right side doesn't exist");
                }
                // check to see what type it is
                if (array[0].equals("seq")) {
                    seqLocation = array[1];
                } else if (array[0].equals("1KG")) {
                    oneKGLocation = array[1];
                } else if (array[0].equals("validation")) {
                    validationDataLocation = array[1];    
                } else {
                    throw new IOException("Line : " + line + ", unknown left side");
                }
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

        @Override
        public void close() {
            // do nothing
        }

        @Override
        public boolean requiresLayout() {
            return false;
        }

        public boolean foundString() {
            return foundString;
        }
    }

}
