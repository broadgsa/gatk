/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.phonehome;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.engine.crypt.CryptUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.io.IOUtils;
import org.broadinstitute.gatk.utils.io.Resource;
import org.broadinstitute.gatk.utils.threading.ThreadEfficiencyMonitor;
import org.jets3t.service.S3Service;
import org.jets3t.service.S3ServiceException;
import org.jets3t.service.impl.rest.httpclient.RestS3Service;
import org.jets3t.service.model.S3Object;
import org.jets3t.service.security.AWSCredentials;
import org.simpleframework.xml.Element;
import org.simpleframework.xml.Serializer;
import org.simpleframework.xml.core.Persister;

import java.io.*;
import java.security.NoSuchAlgorithmException;
import java.security.PublicKey;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;


/**
 * A detailed description of a GATK run, and error if applicable.  Simply create a GATKRunReport
 * with the constructor, providing the walker that was run and the fully instantiated GenomeAnalysisEngine
 * <b>after the run finishes</b> and the GATKRunReport will collect all of the report information
 * into this object.  Call postReport to write out the report, as an XML document, to either STDOUT,
 * a file (in which case the output is gzipped), or with no arguments the report will be posted to the
 * GATK run report database.
 *
 * @author depristo
 * @since 2010
 */
public class GATKRunReport {
    protected static final String REPORT_BUCKET_NAME = "broad.gsa.gatk.run.reports";
    protected static final String TEST_REPORT_BUCKET_NAME = "broad.gsa.gatk.run.reports.test";
    protected final static String AWS_ACCESS_KEY_MD5 = "34d4a26eb2062b3f06e833b28f9a38c6";
    protected final static String AWS_SECRET_KEY_MD5 = "83f2332eec99ef1d7425d5dc5d4b514a";

    private static final DateFormat DATE_FORMAT = new SimpleDateFormat("yyyy/MM/dd HH.mm.ss");

    /**
     * our log
     */
    protected static final Logger logger = Logger.getLogger(GATKRunReport.class);

    /**
     * Default value for the number of milliseconds before an S3 put operation is timed-out.
     * Can be overridden via a constructor argument.
     */
    private static final long S3_DEFAULT_PUT_TIME_OUT_IN_MILLISECONDS = 30 * 1000;

    /**
     * Number of milliseconds before an S3 put operation is timed-out.
     */
    private long s3PutTimeOutInMilliseconds = S3_DEFAULT_PUT_TIME_OUT_IN_MILLISECONDS;

    // -----------------------------------------------------------------
    // elements captured for the report
    // -----------------------------------------------------------------

    @Element(required = false, name = "id")
    private String id;

    @Element(required = false, name = "exception")
    private GATKRunReportException mException;

    @Element(required = true, name = "start-time")
    private String startTime = "ND";

    @Element(required = true, name = "end-time")
    private String endTime;

    @Element(required = true, name = "run-time")
    private long runTime = 0;

    @Element(required = true, name = "walker-name")
    private String walkerName;

    @Element(required = true, name = "svn-version")
    private String svnVersion;

    @Element(required = true, name = "total-memory")
    private long totalMemory;

    @Element(required = true, name = "max-memory")
    private long maxMemory;

    @Element(required = true, name = "user-name")
    private String userName;

    @Element(required = true, name = "host-name")
    private String hostName;

    @Element(required = true, name = "java")
    private String javaVersion;

    @Element(required = true, name = "machine")
    private String machine;

    @Element(required = true, name = "iterations")
    private long nIterations;

    @Element(required = true, name = "tag")
    private String tag;

    @Element(required = true, name = "num-threads")
    private int numThreads;
    @Element(required = true, name = "percent-time-running")
    private String percentTimeRunning;
    @Element(required = true, name = "percent-time-waiting")
    private String percentTimeWaiting;
    @Element(required = true, name = "percent-time-blocking")
    private String percentTimeBlocking;
    @Element(required = true, name = "percent-time-waiting-for-io")
    private String percentTimeWaitingForIO;

    /** The error message, if one occurred, or null if none did */
    public String errorMessage = null;
    /** The error that occurred, if one did, or null if none did */
    public Throwable errorThrown = null;

    /**
     * How should the GATK report its usage?
     */
    public enum PhoneHomeOption {
        /** Disable phone home */
        NO_ET,
        /** Forces the report to go to S3 */
        AWS,
        /** Force output to STDOUT.  For debugging only */
        STDOUT
    }

    /**
     * To allow us to deserial reports from XML
     */
    private GATKRunReport() { }

    /**
     * Read a GATKRunReport from the serialized XML representation in String reportAsXML
     * @param stream an input stream containing a serialized XML report
     * @return a reconstituted GATKRunReport from reportAsXML
     * @throws Exception if parsing fails for any reason
     */
    @Ensures("result != null")
    protected static GATKRunReport deserializeReport(final InputStream stream) throws Exception {
        final Serializer serializer = new Persister();
        return serializer.read(GATKRunReport.class, stream);
    }

    /**
     * Create a new GATKRunReport from a report on S3
     *
     * Assumes that s3Object has already been written to S3, and this function merely
     * fetches it from S3 and deserializes it.  The access keys must have permission to
     * GetObject from S3.
     *
     * @param downloaderAccessKey AWS access key with permission to GetObject from bucketName
     * @param downloaderSecretKey AWS secret key with permission to GetObject from bucketName
     * @param bucketName the name of the bucket holding the report
     * @param s3Object the s3Object we wrote to S3 in bucketName that we want to get back and decode
     * @return a deserialized report derived from s3://bucketName/s3Object.getName()
     * @throws Exception
     */
    @Ensures("result != null")
    protected static GATKRunReport deserializeReport(final String downloaderAccessKey,
                                                     final String downloaderSecretKey,
                                                     final String bucketName,
                                                     final S3Object s3Object) throws Exception {
        final S3Service s3Service = initializeAWSService(downloaderAccessKey, downloaderSecretKey);

        // Retrieve the whole data object we created previously
        final S3Object objectComplete = s3Service.getObject(bucketName, s3Object.getName());

        // Read the data from the object's DataInputStream using a loop, and print it out.
        return deserializeReport(new GZIPInputStream(objectComplete.getDataInputStream()));
    }

    /**
     * Create a new RunReport and population all of the fields with values from the walker and engine.
     * Allows the S3 put timeout to be explicitly set.
     *
     * @param walker the GATK walker that we ran
     * @param e the exception caused by running this walker, or null if we completed successfully
     * @param engine the GAE we used to run the walker, so we can fetch runtime, args, etc
     * @param type the GATK phone home setting
     * @param s3PutTimeOutInMilliseconds number of milliseconds to wait before timing out an S3 put operation
     */
    public GATKRunReport(final Walker<?,?> walker, final Exception e, final GenomeAnalysisEngine engine, final PhoneHomeOption type,
                         final long s3PutTimeOutInMilliseconds) {
        this(walker, e, engine, type);
        this.s3PutTimeOutInMilliseconds = s3PutTimeOutInMilliseconds;
    }

    /**
     * Create a new RunReport and population all of the fields with values from the walker and engine.
     * Leaves the S3 put timeout set to the default value of S3_DEFAULT_PUT_TIME_OUT_IN_MILLISECONDS.
     *
     * @param walker the GATK walker that we ran
     * @param e the exception caused by running this walker, or null if we completed successfully
     * @param engine the GAE we used to run the walker, so we can fetch runtime, args, etc
     * @param type the GATK phone home setting
     */
    public GATKRunReport(final Walker<?,?> walker, final Exception e, final GenomeAnalysisEngine engine, final PhoneHomeOption type) {
        if ( type == PhoneHomeOption.NO_ET )
            throw new ReviewedGATKException("Trying to create a run report when type is NO_ET!");

        logger.debug("Aggregating data for run report");

        // what did we run?
        id = org.apache.commons.lang.RandomStringUtils.randomAlphanumeric(32);
        walkerName = engine.getWalkerName(walker.getClass());
        svnVersion = CommandLineGATK.getVersionNumber();

        // runtime performance metrics
        Date end = new java.util.Date();
        endTime = DATE_FORMAT.format(end);
        if ( engine.getStartTime() != null ) { // made it this far during initialization
            startTime = DATE_FORMAT.format(engine.getStartTime());
            runTime = (end.getTime() - engine.getStartTime().getTime()) / 1000L; // difference in seconds
        }

        // deal with memory usage
        Runtime.getRuntime().gc(); // call GC so totalMemory is ~ used memory
        maxMemory = Runtime.getRuntime().maxMemory();
        totalMemory = Runtime.getRuntime().totalMemory();

        // we can only do some operations if an error hasn't occurred
        if ( engine.getCumulativeMetrics() != null ) {
            // it's possible we aborted so early that these data structures arent initialized
            nIterations = engine.getCumulativeMetrics().getNumIterations();
        }

        tag = engine.getArguments().tag;

        // user and hostname -- information about the runner of the GATK
        userName = System.getProperty("user.name");
        hostName = Utils.resolveHostname();

        // basic java information
        javaVersion = Utils.join("-", Arrays.asList(System.getProperty("java.vendor"), System.getProperty("java.version")));
        machine = Utils.join("-", Arrays.asList(System.getProperty("os.name"), System.getProperty("os.arch")));

        // if there was an exception, capture it
        this.mException = e == null ? null : new GATKRunReportException(e);

        numThreads = engine.getTotalNumberOfThreads();
        percentTimeRunning = getThreadEfficiencyPercent(engine, ThreadEfficiencyMonitor.State.USER_CPU);
        percentTimeBlocking = getThreadEfficiencyPercent(engine, ThreadEfficiencyMonitor.State.BLOCKING);
        percentTimeWaiting = getThreadEfficiencyPercent(engine, ThreadEfficiencyMonitor.State.WAITING);
        percentTimeWaitingForIO = getThreadEfficiencyPercent(engine, ThreadEfficiencyMonitor.State.WAITING_FOR_IO);
    }

    /**
     * Get the random alpha-numeric ID of this GATKRunReport
     * @return a non-null string ID
     */
    @Ensures("result != null")
    public String getID() {
        return id;
    }

    /**
     * Return a string representing the percent of time the GATK spent in state, if possible.  Otherwise return NA
     *
     * @param engine the GATK engine whose threading efficiency info we will use
     * @param state the state whose occupancy we wish to know
     * @return a string representation of the percent occupancy of state, or NA is not possible
     */
    @Requires({"engine != null", "state != null"})
    @Ensures("result != null")
    private String getThreadEfficiencyPercent(final GenomeAnalysisEngine engine, final ThreadEfficiencyMonitor.State state) {
        final ThreadEfficiencyMonitor tem = engine.getThreadEfficiencyMonitor();
        return tem == null ? "NA" : String.format("%.2f", tem.getStatePercent(state));
    }

    /**
     * Get a filename (no path) appropriate for this report
     *
     * @return a non-null string filename
     */
    @Ensures("result != null")
    protected String getReportFileName() {
        return getID() + ".report.xml.gz";
    }

    // ---------------------------------------------------------------------------
    //
    // Main public interface method for posting reports
    //
    // ---------------------------------------------------------------------------

    /**
     * Post this GATK report to the destination implied by the PhoneHomeOption type
     *
     * Guaranteed to never throw an exception (exception noted below) and to return
     * with a reasonable (~10 seconds) time regardless of successful writing of the report.
     *
     * @throws IllegalArgumentException if type == null
     * @param type the type of phoning home we want to do
     * @return true if a report was successfully written, false otherwise
     */
    public boolean postReport(final PhoneHomeOption type) {
        if ( type == null ) throw new IllegalArgumentException("type cannot be null");

        logger.debug("Posting report of type " + type);
        switch (type) {
            case NO_ET: // don't do anything
                return false;
            case AWS:
                wentToAWS = true;
                return postReportToAWSS3() != null;
            case STDOUT:
                return postReportToStream(System.out);
            default:
                exceptDuringRunReport("BUG: unexpected PhoneHomeOption ");
                return false;
        }
    }

    // ---------------------------------------------------------------------------
    //
    // Code for sending reports to local files
    //
    // ---------------------------------------------------------------------------

    /**
     * Write an XML representation of this report to the stream, throwing a GATKException if the marshalling
     * fails for any reason.
     *
     * @param stream an output stream to write the report to
     */
    @Requires("stream != null")
    protected boolean postReportToStream(final OutputStream stream) {
        final Serializer serializer = new Persister();
        try {
            serializer.write(this, stream);
            return true;
        } catch (Exception e) {
            return false;
        }
    }

    // ---------------------------------------------------------------------------
    //
    // Code for sending reports to s3
    //
    // ---------------------------------------------------------------------------

    /**
     * Get the name of the S3 bucket where we should upload this report
     *
     * @return the string name of the s3 bucket
     */
    @Ensures("result != null")
    protected String getS3ReportBucket() {
        return s3ReportBucket;
    }

    /**
     * Decrypts encrypted AWS key from encryptedKeySource
     * @param encryptedKeySource a file containing an encrypted AWS key
     * @return a decrypted AWS key as a String
     */
    @Ensures("result != null")
    public static String decryptAWSKey(final File encryptedKeySource) throws FileNotFoundException {
        if ( encryptedKeySource == null ) throw new IllegalArgumentException("encryptedKeySource cannot be null");
        return decryptAWSKey(new FileInputStream(encryptedKeySource));
    }

    /**
     * @see #decryptAWSKey(java.io.File) but with input from an inputstream
     */
    @Requires("encryptedKeySource != null")
    @Ensures("result != null")
    private static String decryptAWSKey(final InputStream encryptedKeySource) {
        final PublicKey key = CryptUtils.loadGATKDistributedPublicKey();
        final byte[] fromDisk = IOUtils.readStreamIntoByteArray(encryptedKeySource);
        final byte[] decrypted = CryptUtils.decryptData(fromDisk, key);
        return new String(decrypted);
    }

    /**
     * Get the decrypted AWS key sorted in the resource directories of name
     * @param name the name of the file containing the needed AWS key
     * @return a non-null GATK
     */
    @Requires("name != null")
    @Ensures("result != null")
    private static String getAWSKey(final String name) {
        final Resource resource = new Resource(name, GATKRunReport.class);
        return decryptAWSKey(resource.getResourceContentsAsStream());
    }

    /**
     * Get the AWS access key for the GATK user
     * @return a non-null AWS access key for the GATK user
     */
    @Ensures("result != null")
    protected static String getAWSUploadAccessKey() {
        return getAWSKey("resources/GATK_AWS_access.key");
    }

    /**
     * Get the AWS secret key for the GATK user
     * @return a non-null AWS secret key for the GATK user
     */
    @Ensures("result != null")
    protected static String getAWSUploadSecretKey() {
        return getAWSKey("resources/GATK_AWS_secret.key");
    }

    /**
     * Check that the AWS keys can be decrypted and are what we expect them to be
     *
     * @throws ReviewedGATKException if anything goes wrong
     */
    public static void checkAWSAreValid() {
        try {
            final String accessKeyMD5 = Utils.calcMD5(getAWSUploadAccessKey());
            final String secretKeyMD5 = Utils.calcMD5(getAWSUploadSecretKey());

            if ( ! AWS_ACCESS_KEY_MD5.equals(accessKeyMD5) ) {
                throw new ReviewedGATKException("Invalid AWS access key found, expected MD5 " + AWS_ACCESS_KEY_MD5 + " but got " + accessKeyMD5);
            }
            if ( ! AWS_SECRET_KEY_MD5.equals(secretKeyMD5) ) {
                throw new ReviewedGATKException("Invalid AWS secret key found, expected MD5 " + AWS_SECRET_KEY_MD5 + " but got " + secretKeyMD5);
            }

        } catch ( Exception e ) {
            throw new ReviewedGATKException("Couldn't decrypt AWS keys, something is wrong with the GATK distribution");
        }
    }

    /**
     * Get an initialized S3Service for use in communicating with AWS/s3
     *
     * @param awsAccessKey our AWS access key to use
     * @param awsSecretKey our AWS secret key to use
     * @return an initialized S3Service object that can be immediately used to interact with S3
     * @throws S3ServiceException
     */
    @Requires({"awsAccessKey != null", "awsSecretKey != null"})
    @Ensures("result != null")
    protected static S3Service initializeAWSService(final String awsAccessKey, final String awsSecretKey) throws S3ServiceException {
        // To communicate with S3, create a class that implements an S3Service. We will use the REST/HTTP
        // implementation based on HttpClient, as this is the most robust implementation provided with JetS3t.
        final AWSCredentials awsCredentials = new AWSCredentials(awsAccessKey, awsSecretKey);
        return new RestS3Service(awsCredentials);
    }

    /**
     * A runnable that pushes this GATKReport up to s3.
     *
     * Should be run in a separate thread so we can time it out if something is taking too long
     */
    private class S3PutRunnable implements Runnable {
        /** Was the upload operation successful? */
        public final AtomicBoolean isSuccess;
        /** The name of this report */
        private final String filename;
        /** The contents of this report */
        private final byte[] contents;

        /** The s3Object that we created to upload, or null if it failed */
        public S3Object s3Object = null;

        @Requires({"filename != null", "contents != null"})
        public S3PutRunnable(final String filename, final byte[] contents){
            this.isSuccess = new AtomicBoolean();
            this.filename = filename;
            this.contents = contents;
        }

        public void run() {
            try {
                switch ( awsMode ) {
                    case FAIL_WITH_EXCEPTION:
                        throw new IllegalStateException("We are throwing an exception for testing purposes");
                    case TIMEOUT:
                        try {
                            Thread.sleep(s3PutTimeOutInMilliseconds * 100);
                        } catch ( InterruptedException e ) {
                            // supposed to be empty
                        }
                        break;
                    case NORMAL:
                        // IAM GATK user credentials -- only right is to PutObject into broad.gsa.gatk.run.reports bucket
                        final S3Service s3Service = initializeAWSService(getAWSUploadAccessKey(), getAWSUploadSecretKey());

                        // Create an S3Object based on a file, with Content-Length set automatically and
                        // Content-Type set based on the file's extension (using the Mimetypes utility class)
                        final S3Object fileObject = new S3Object(filename, contents);
                        //logger.info("Created S3Object" + fileObject);
                        //logger.info("Uploading " + localFile + " to AWS bucket");
                        s3Object = s3Service.putObject(getS3ReportBucket(), fileObject);
                        isSuccess.set(true);
                        break;
                    default:
                        throw new IllegalStateException("Unexpected AWS exception");
                }
            } catch ( S3ServiceException e ) {
                exceptDuringRunReport("S3 exception occurred", e);
            } catch ( NoSuchAlgorithmException e ) {
                exceptDuringRunReport("Couldn't calculate MD5", e);
            } catch ( IOException e ) {
                exceptDuringRunReport("Couldn't read report file", e);
            } catch ( Exception e ) {
                exceptDuringRunReport("An unexpected exception occurred during posting", e);
            }
        }
    }

    /**
     * Post this GATK report to the AWS s3 GATK_Run_Report log
     *
     * @return the s3Object pointing to our pushed report, or null if we failed to push
     */
    protected S3Object postReportToAWSS3() {
        // modifying example code from http://jets3t.s3.amazonaws.com/toolkit/code-samples.html
        this.hostName = Utils.resolveHostname(); // we want to fill in the host name
        final String key = getReportFileName();
        logger.debug("Generating GATK report to AWS S3 with key " + key);

        try {
            // create an byte output stream so we can capture the output as a byte[]
            final ByteArrayOutputStream byteStream = new ByteArrayOutputStream(8096);
            final OutputStream outputStream = new GZIPOutputStream(byteStream);
            postReportToStream(outputStream);
            outputStream.close();
            final byte[] report = byteStream.toByteArray();

            // stop us from printing the annoying, and meaningless, mime types warning
            final Logger mimeTypeLogger = Logger.getLogger(org.jets3t.service.utils.Mimetypes.class);
            mimeTypeLogger.setLevel(Level.FATAL);

            // Set the S3 upload on its own thread with timeout:
            final S3PutRunnable s3run = new S3PutRunnable(key,report);
            final Thread s3thread = new Thread(s3run);
            s3thread.setDaemon(true);
            s3thread.setName("S3Put-Thread");
            s3thread.start();

            s3thread.join(s3PutTimeOutInMilliseconds);

            if(s3thread.isAlive()){
                s3thread.interrupt();
                exceptDuringRunReport("Run statistics report upload to AWS S3 timed-out");
            } else if(s3run.isSuccess.get()) {
                logger.info("Uploaded run statistics report to AWS S3");
                logger.debug("Uploaded to AWS: " + s3run.s3Object);
                return s3run.s3Object;
            } else {
                // an exception occurred, the thread should have already invoked the exceptDuringRunReport function
            }
        } catch ( IOException e ) {
            exceptDuringRunReport("Couldn't read report file", e);
        } catch ( InterruptedException e) {
            exceptDuringRunReport("Run statistics report upload interrupted", e);
        }

        return null;
    }

    // ---------------------------------------------------------------------------
    //
    // Error handling code
    //
    // ---------------------------------------------------------------------------

    /**
     * Note that an exception occurred during creating or writing this report
     * @param msg the message to print
     * @param e the exception that occurred
     */
    @Ensures("exceptionOccurredDuringPost()")
    private void exceptDuringRunReport(final String msg, final Throwable e) {
        this.errorMessage = msg;
        this.errorThrown = e;
        logger.debug("A problem occurred during GATK run reporting [*** everything is fine, but no report could be generated; please do not post this to the support forum ***].  Message is: " + msg + ".  Error message is: " + e.getMessage());
    }

    /**
     * Note that an exception occurred during creating or writing this report
     * @param msg the message to print
     */
    @Ensures("exceptionOccurredDuringPost()")
    private void exceptDuringRunReport(final String msg) {
        this.errorMessage = msg;
        logger.debug("A problem occurred during GATK run reporting [*** everything is fine, but no report could be generated; please do not post this to the support forum ***].  Message is " + msg);
    }

    /**
     * Did an error occur during the posting of this run report?
     * @return true if so, false if not
     */
    public boolean exceptionOccurredDuringPost() {
        return getErrorMessage() != null;
    }

    /**
     * If an error occurred during posting of this report, retrieve the message of the error that occurred, or null if
     * no error occurred
     * @return a string describing the error that occurred, or null if none did
     */
    public String getErrorMessage() {
        return errorMessage;
    }

    /**
     * Get the throwable that caused the exception during posting of this message, or null if none was available
     *
     * Note that getting a null valuable from this function doesn't not imply that no error occurred.  Some
     * errors that occurred many not have generated a throwable.
     *
     * @return the Throwable that caused the error, or null if no error occurred or was not caused by a throwable
     */
    public Throwable getErrorThrown() {
        return errorThrown;
    }

    /**
     * Helper method to format the exception that occurred during posting, or a string saying none occurred
     * @return a non-null string
     */
    @Ensures("result != null")
    protected String formatError() {
        return exceptionOccurredDuringPost()
                ? String.format("Exception message=%s with cause=%s", getErrorMessage(), getErrorThrown())
                : "No exception occurred";
    }

    // ---------------------------------------------------------------------------
    //
    // Equals and hashcode -- purely for comparing reports for testing
    //
    // ---------------------------------------------------------------------------

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        GATKRunReport that = (GATKRunReport) o;

        if (maxMemory != that.maxMemory) return false;
        if (nIterations != that.nIterations) return false;
        if (numThreads != that.numThreads) return false;
        if (runTime != that.runTime) return false;
        if (totalMemory != that.totalMemory) return false;
        if (endTime != null ? !endTime.equals(that.endTime) : that.endTime != null) return false;
        if (hostName != null ? !hostName.equals(that.hostName) : that.hostName != null) return false;
        if (id != null ? !id.equals(that.id) : that.id != null) return false;
        if (javaVersion != null ? !javaVersion.equals(that.javaVersion) : that.javaVersion != null) return false;
        if (mException != null ? !mException.equals(that.mException) : that.mException != null) return false;
        if (machine != null ? !machine.equals(that.machine) : that.machine != null) return false;
        if (percentTimeBlocking != null ? !percentTimeBlocking.equals(that.percentTimeBlocking) : that.percentTimeBlocking != null)
            return false;
        if (percentTimeRunning != null ? !percentTimeRunning.equals(that.percentTimeRunning) : that.percentTimeRunning != null)
            return false;
        if (percentTimeWaiting != null ? !percentTimeWaiting.equals(that.percentTimeWaiting) : that.percentTimeWaiting != null)
            return false;
        if (percentTimeWaitingForIO != null ? !percentTimeWaitingForIO.equals(that.percentTimeWaitingForIO) : that.percentTimeWaitingForIO != null)
            return false;
        if (startTime != null ? !startTime.equals(that.startTime) : that.startTime != null) return false;
        if (svnVersion != null ? !svnVersion.equals(that.svnVersion) : that.svnVersion != null) return false;
        if (tag != null ? !tag.equals(that.tag) : that.tag != null) return false;
        if (userName != null ? !userName.equals(that.userName) : that.userName != null) return false;
        if (walkerName != null ? !walkerName.equals(that.walkerName) : that.walkerName != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = id != null ? id.hashCode() : 0;
        result = 31 * result + (mException != null ? mException.hashCode() : 0);
        result = 31 * result + (startTime != null ? startTime.hashCode() : 0);
        result = 31 * result + (endTime != null ? endTime.hashCode() : 0);
        result = 31 * result + (int) (runTime ^ (runTime >>> 32));
        result = 31 * result + (walkerName != null ? walkerName.hashCode() : 0);
        result = 31 * result + (svnVersion != null ? svnVersion.hashCode() : 0);
        result = 31 * result + (int) (totalMemory ^ (totalMemory >>> 32));
        result = 31 * result + (int) (maxMemory ^ (maxMemory >>> 32));
        result = 31 * result + (userName != null ? userName.hashCode() : 0);
        result = 31 * result + (hostName != null ? hostName.hashCode() : 0);
        result = 31 * result + (javaVersion != null ? javaVersion.hashCode() : 0);
        result = 31 * result + (machine != null ? machine.hashCode() : 0);
        result = 31 * result + (int) (nIterations ^ (nIterations >>> 32));
        result = 31 * result + (tag != null ? tag.hashCode() : 0);
        result = 31 * result + numThreads;
        result = 31 * result + (percentTimeRunning != null ? percentTimeRunning.hashCode() : 0);
        result = 31 * result + (percentTimeWaiting != null ? percentTimeWaiting.hashCode() : 0);
        result = 31 * result + (percentTimeBlocking != null ? percentTimeBlocking.hashCode() : 0);
        result = 31 * result + (percentTimeWaitingForIO != null ? percentTimeWaitingForIO.hashCode() : 0);
        return result;
    }

    // ---------------------------------------------------------------------------
    //
    // Code specifically for testing the GATKRunReport
    //
    // ---------------------------------------------------------------------------

    /**
     * Enum specifying how the S3 uploader should behave.  Must be normal by default.  Purely for testing purposes
     */
    protected enum AWSMode {
        NORMAL, // write normally to AWS
        FAIL_WITH_EXCEPTION, // artificially fail during writing
        TIMEOUT // sleep, so we time out
    }
    /** Our AWS mode */
    private AWSMode awsMode = AWSMode.NORMAL;
    /** The bucket were we send the GATK report on AWS/s3 */
    private String s3ReportBucket = REPORT_BUCKET_NAME;
    /** Did we send the report to AWS? */
    private boolean wentToAWS = false;

    /**
     * Send the report to the AWS test bucket -- for testing only
     */
    protected void sendAWSToTestBucket() {
        s3ReportBucket = TEST_REPORT_BUCKET_NAME;
    }

    /**
     * Has the report been written to AWS?
     *
     * Does not imply anything about the success of the send, just that it was attempted
     *
     * @return true if the report has been sent to AWS, false otherwise
     */
    protected boolean wentToAWS() {
        return wentToAWS;
    }

    /**
     * Purely for testing purposes.  Tells the AWS uploader whether to actually upload or simulate errors
     * @param mode what we want to do
     */
    @Requires("mode != null")
    protected void setAwsMode(final AWSMode mode) {
        this.awsMode = mode;
    }
}
