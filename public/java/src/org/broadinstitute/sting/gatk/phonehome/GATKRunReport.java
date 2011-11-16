/*
 * Copyright (c) 2010, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.phonehome;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.jets3t.service.S3Service;
import org.jets3t.service.S3ServiceException;
import org.jets3t.service.impl.rest.httpclient.RestS3Service;
import org.jets3t.service.model.S3Object;
import org.jets3t.service.security.AWSCredentials;
import org.simpleframework.xml.Element;
import org.simpleframework.xml.ElementList;
import org.simpleframework.xml.Serializer;
import org.simpleframework.xml.core.Persister;
import org.simpleframework.xml.stream.Format;
import org.simpleframework.xml.stream.HyphenStyle;

import java.io.*;
import java.security.NoSuchAlgorithmException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.zip.GZIPOutputStream;


/**
 * @author depristo
 *
 * A detailed description of a GATK run, and error if applicable.  Simply create a GATKRunReport
 * with the constructor, providing the walker that was run and the fully instantiated GenomeAnalysisEngine
 * <b>after the run finishes</b> and the GATKRunReport will collect all of the report information
 * into this object.  Call postReport to write out the report, as an XML document, to either STDOUT,
 * a file (in which case the output is gzipped), or with no arguments the report will be posted to the
 * GATK run report database.
 */
public class GATKRunReport {
    /**
     * The root file system directory where we keep common report data
     */
    private static File REPORT_DIR = new File("/humgen/gsa-hpprojects/GATK/reports");

    private static final String REPORT_BUCKET_NAME = "GATK_Run_Reports";

    /**
     * The full path to the direct where submitted (and uncharacterized) report files are written
     */
    private static File REPORT_SUBMIT_DIR = new File(REPORT_DIR.getAbsolutePath() + "/submitted");

    /**
     * Full path to the sentinel file that controls whether reports are written out.  If this file doesn't
     * exist, no long will be written
     */
    private static File REPORT_SENTINEL = new File(REPORT_DIR.getAbsolutePath() + "/ENABLE");

    /**
     * our log
     */
    protected static Logger logger = Logger.getLogger(GATKRunReport.class);


    // the listing of the fields is somewhat important; this is the order that the simple XML will output them
    @ElementList(required = true, name = "gatk_header_Information")
    private List<String> mGATKHeader;

    @Element(required = false, name = "id")
    private final String id;

    @Element(required = false, name = "exception")
    private final ExceptionToXML mException;

    @Element(required = true, name = "working_directory")
    private String currentPath;

    @Element(required = true, name = "start_time")
    private String startTime = "ND";

    @Element(required = true, name = "end_time")
    private String endTime;

    @Element(required = true, name = "run_time")
    private long runTime = 0;

    @Element(required = true, name = "command_line")
    private String cmdLine = "COULD NOT BE DETERMINED";

    @Element(required = true, name = "walker_name")
    private String walkerName;

    @Element(required = true, name = "svn_version")
    private String svnVersion;

    @Element(required = true, name = "total_memory")
    private long totalMemory;

    @Element(required = true, name = "max_memory")
    private long maxMemory;

    @Element(required = true, name = "java_tmp_directory")
    private String tmpDir;

    @Element(required = true, name = "user_name")
    private String userName;

    @Element(required = true, name = "host_name")
    private String hostName;

    @Element(required = true, name = "java")
    private String java;

    @Element(required = true, name = "machine")
    private String machine;

    @Element(required = true, name = "iterations")
    private long nIterations;

    @Element(required = true, name = "reads")
    private long nReads;

    public enum PhoneHomeOption {
        /** Disable phone home */
        NO_ET,
        /** Standard option.  Writes to local repository if it can be found, or S3 otherwise */
        STANDARD,
        /** Force output to STDOUT.  For debugging only */
        STDOUT,
        /** Force output to S3.  For debugging only */
        AWS_S3   // todo -- remove me -- really just for testing purposes
    }

    private static final DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH.mm.ss");

    /**
     * Create a new RunReport and population all of the fields with values from the walker and engine
     *
     * @param walker the GATK walker that we ran
     * @param e the exception caused by running this walker, or null if we completed successfully
     * @param engine the GAE we used to run the walker, so we can fetch runtime, args, etc
     */
    public GATKRunReport(Walker<?,?> walker, Exception e, GenomeAnalysisEngine engine, PhoneHomeOption type) {
        if ( type == PhoneHomeOption.NO_ET )
            throw new ReviewedStingException("Trying to create a run report when type is NO_ET!");

        logger.debug("Aggregating data for run report");

        mGATKHeader = CommandLineGATK.createApplicationHeader();
        currentPath = System.getProperty("user.dir");

        // what did we run?
        id = org.apache.commons.lang.RandomStringUtils.randomAlphanumeric(32);
        try {
            cmdLine = engine.createApproximateCommandLineArgumentString(engine, walker);
        } catch (Exception ignore) { }

        walkerName = engine.getWalkerName(walker.getClass());
        svnVersion = CommandLineGATK.getVersionNumber();

        // runtime performance metrics
        Date end = new java.util.Date();
        endTime = dateFormat.format(end);
        if ( engine.getStartTime() != null ) { // made it this far during initialization
            startTime = dateFormat.format(engine.getStartTime());
            runTime = (end.getTime() - engine.getStartTime().getTime()) / 1000L; // difference in seconds
        }
        tmpDir = System.getProperty("java.io.tmpdir");

        // deal with memory usage
        Runtime.getRuntime().gc(); // call GC so totalMemory is ~ used memory
        maxMemory = Runtime.getRuntime().maxMemory();
        totalMemory = Runtime.getRuntime().totalMemory();

        // we can only do some operations if an error hasn't occurred
        if ( engine.getCumulativeMetrics() != null ) {
            // it's possible we aborted so early that these data structures arent initialized
            nIterations = engine.getCumulativeMetrics().getNumIterations();
            nReads = engine.getCumulativeMetrics().getNumReadsSeen();
        }

        // user and hostname -- information about the runner of the GATK
        userName = System.getProperty("user.name");
        hostName = "unknown"; // resolveHostname();

        // basic java information
        java = Utils.join("-", Arrays.asList(System.getProperty("java.vendor"), System.getProperty("java.version")));
        machine = Utils.join("-", Arrays.asList(System.getProperty("os.name"), System.getProperty("os.arch")));

        // if there was an exception, capture it
        this.mException = e == null ? null : new ExceptionToXML(e);
    }

    public String getID() {
        return id;
    }


    public void postReport(PhoneHomeOption type) {
        logger.debug("Posting report of type " + type);
        switch (type) {
            case NO_ET: // don't do anything
                break;
            case STANDARD:
                if ( repositoryIsOnline() ) {
                    postReportToLocalDisk(REPORT_SUBMIT_DIR);
                } else {
                    postReportToAWSS3();
                }
                break;
            case STDOUT:
                postReportToStream(System.out);
                break;
            case AWS_S3:
                postReportToAWSS3();
                break;
            default:
                exceptDuringRunReport("BUG: unexcepted PhoneHomeOption ");
                break;
        }
    }

    /**
     * Write an XML representation of this report to the stream, throwing a StingException if the marshalling
     * fails for any reason.
     *
     * @param stream
     */
    private void postReportToStream(OutputStream stream) {
        Serializer serializer = new Persister(new Format(new HyphenStyle()));
        try {
            serializer.write(this, stream);
            //throw new StingException("test");
        } catch (Exception e) {
            throw new ReviewedStingException("Failed to marshal the data to the file " + stream, e);
        }
    }

    /**
     * Opens the destination file and writes a gzipped version of the XML report there.
     *
     * @param destination
     * @throws IOException
     */
    private void postReportToFile(File destination) throws IOException {
        BufferedOutputStream out =
                new BufferedOutputStream(
                        new GZIPOutputStream(
                                new FileOutputStream(destination)));
        try {
            postReportToStream(out);
        } finally {
            out.close();
        }
    }

    /**
     * Main entry point to writing reports to disk.  Posts the XML report to the common GATK run report repository.
     * If this process fails for any reason, all exceptions are handled and this routine merely prints a warning.
     * That is, postReport() is guarenteed not to fail for any reason.
     */
    private File postReportToLocalDisk(File rootDir) {
        String filename = getID() + ".report.xml.gz";
        File file = new File(rootDir, filename);
        try {
            postReportToFile(file);
            logger.debug("Wrote report to " + file);
            return file;
        } catch ( Exception e ) {
            // we catch everything, and no matter what eat the error
            exceptDuringRunReport("Couldn't read report file", e);
            file.delete();
            return null;
        }
    }

    private void postReportToAWSS3() {
        // modifying example code from http://jets3t.s3.amazonaws.com/toolkit/code-samples.html
        this.hostName = Utils.resolveHostname(); // we want to fill in the host name
        File localFile = postReportToLocalDisk(new File("./"));
        logger.debug("Generating GATK report to AWS S3 based on local file " + localFile);
        if ( localFile != null ) { // we succeeded in creating the local file
            localFile.deleteOnExit();
            try {
                // stop us from printing the annoying, and meaningless, mime types warning
                Logger mimeTypeLogger = Logger.getLogger(org.jets3t.service.utils.Mimetypes.class);
                mimeTypeLogger.setLevel(Level.FATAL);

                // Your Amazon Web Services (AWS) login credentials are required to manage S3 accounts. These credentials
                // are stored in an AWSCredentials object:

                // IAM GATK user credentials -- only right is to PutObject into GATK_Run_Report bucket
                String awsAccessKey = "AKIAJXU7VIHBPDW4TDSQ"; // GATK AWS user
                String awsSecretKey = "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"; // GATK AWS user
                AWSCredentials awsCredentials = new AWSCredentials(awsAccessKey, awsSecretKey);

                // To communicate with S3, create a class that implements an S3Service. We will use the REST/HTTP
                // implementation based on HttpClient, as this is the most robust implementation provided with JetS3t.
                S3Service s3Service = new RestS3Service(awsCredentials);

                // Create an S3Object based on a file, with Content-Length set automatically and
                // Content-Type set based on the file's extension (using the Mimetypes utility class)
                S3Object fileObject = new S3Object(localFile);
                //logger.info("Created S3Object" + fileObject);
                //logger.info("Uploading " + localFile + " to AWS bucket");
                S3Object s3Object = s3Service.putObject(REPORT_BUCKET_NAME, fileObject);
                logger.debug("Uploaded to AWS: " + s3Object);
                logger.info("Uploaded run statistics report to AWS S3");
            } catch ( S3ServiceException e ) {
                exceptDuringRunReport("S3 exception occurred", e);
            } catch ( NoSuchAlgorithmException e ) {
                exceptDuringRunReport("Couldn't calculate MD5", e);
            } catch ( IOException e ) {
                exceptDuringRunReport("Couldn't read report file", e);
            }
        }
    }

    private void exceptDuringRunReport(String msg, Throwable e) {
        logger.debug("A problem occurred during GATK run reporting [*** everything is fine, but no report could be generated; please do not post this to the support forum ***].  Message is: " + msg + ".  Error message is: " + e.getMessage());
        //e.printStackTrace();
    }

    private void exceptDuringRunReport(String msg) {
        logger.debug("A problem occurred during GATK run reporting [*** everything is fine, but no report could be generated; please do not post this to the support forum ***].  Message is " + msg);
    }


    /**
     * Returns true if and only if the common run report repository is available and online to receive reports
     *
     * @return
     */
    private boolean repositoryIsOnline() {
        return REPORT_SENTINEL.exists();
    }

    /**
     * A helper class for formatting in XML the throwable chain starting at e.
     */
    private class ExceptionToXML {
        @Element(required = false, name = "message")
        String message = null;

        @ElementList(required = false, name = "stacktrace")
        final List<String> stackTrace = new ArrayList<String>();

        @Element(required = false, name = "cause")
        ExceptionToXML cause = null;

        @Element(required = false, name = "is-user-exception")
        Boolean isUserException;

        public ExceptionToXML(Throwable e) {
            message = e.getMessage();
            isUserException = e instanceof UserException;
            for (StackTraceElement element : e.getStackTrace()) {
                stackTrace.add(element.toString());
            }

            if ( e.getCause() != null ) {
                cause = new ExceptionToXML(e.getCause());
            }
        }
    }
}
