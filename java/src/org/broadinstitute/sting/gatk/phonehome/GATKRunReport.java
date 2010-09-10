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

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.CommandLineUtils;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GATKException;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.simpleframework.xml.Element;
import org.simpleframework.xml.ElementList;
import org.simpleframework.xml.Serializer;
import org.simpleframework.xml.core.Persister;
import org.simpleframework.xml.stream.Format;
import org.simpleframework.xml.stream.HyphenStyle;

import java.io.*;
import java.net.InetAddress;
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

    @Element(required = true, name = "argument_collection")
    private final GATKArgumentCollection mCollection;

    @Element(required = true, name = "working_directory")
    private String currentPath;

    @Element(required = true, name = "start_time")
    private String startTime;

    @Element(required = true, name = "end_time")
    private String endTime;

    @Element(required = true, name = "run_time")
    private long runTime;

    @Element(required = true, name = "command_line")
    private String cmdLine;

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

//    @Element(required = true, name = "read_metrics")
//    private String readMetrics;

    // TODO
    // todo md5 all filenames
    // todo size of filenames

    public enum PhoneHomeOption {
        NO_ET,
        STANDARD,
        DEV,
        STDOUT
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
            throw new GATKException("Trying to create a run report when type is NO_ET!");

        mGATKHeader = CommandLineGATK.createApplicationHeader();
        currentPath = System.getProperty("user.dir");

        // what did we run?
        id = org.apache.commons.lang.RandomStringUtils.randomAlphanumeric(32);
        cmdLine = CommandLineUtils.createApproximateCommandLineArgumentString(engine, walker);
        this.mCollection = engine.getArguments();
        walkerName = engine.getWalkerName(walker.getClass());
        svnVersion = CommandLineGATK.getVersionNumber();

        // runtime performance metrics
        startTime = dateFormat.format(engine.getStartTime());
        Date end = new java.util.Date();
        endTime = dateFormat.format(end);
        runTime = (end.getTime() - engine.getStartTime().getTime()) / 1000L; // difference in seconds
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
        hostName = resolveHostname();

        // basic java information
        java = Utils.join("-", Arrays.asList(System.getProperty("java.vendor"), System.getProperty("java.version")));
        machine = Utils.join("-", Arrays.asList(System.getProperty("os.name"), System.getProperty("os.arch")));

        // if there was an exception, capture it
        this.mException = e == null ? null : new ExceptionToXML(e);
    }

    public String getID() {
        return id;
    }


    /**
     * Helper utility that calls into the InetAddress system to resolve the hostname.  If this fails,
     * unresolvable gets returned instead.
     *
     * @return
     */
    private String resolveHostname() {
        try {
            return InetAddress.getLocalHost().getCanonicalHostName();
        }
        catch (java.net.UnknownHostException uhe) { // [beware typo in code sample -dmw]
            return "unresolvable";
            // handle exception
        }
    }

    /**
     * Write an XML representation of this report to the stream, throwing a StingException if the marshalling
     * fails for any reason.
     *
     * @param stream
     */
    public void postReport(OutputStream stream) {
        Serializer serializer = new Persister(new Format(new HyphenStyle()));
        try {
            serializer.write(this, stream);
            //throw new StingException("test");
        } catch (Exception e) {
            throw new GATKException("Failed to marshal the data to the file " + stream, e);
        }
    }

    /**
     * Opens the destination file and writes a gzipped version of the XML report there.
     *
     * @param destination
     * @throws IOException
     */
    public void postReport(File destination) throws IOException {
        BufferedOutputStream out =
                new BufferedOutputStream(
                        new GZIPOutputStream(
                                new FileOutputStream(destination)));
        try {
            postReport(out);
        } finally {
            out.close();
        }
    }

    /**
     * Main entry point to writing reports to disk.  Posts the XML report to the common GATK run report repository.
     * If this process fails for any reason, all exceptions are handled and this routine merely prints a warning.
     * That is, postReport() is guarenteed not to fail for any reason.
     */
    public void postReport() {
        try {
            if ( repositoryIsOnline() ) {
                String filename = getID() + ".report.xml.gz";
                File file = new File(REPORT_SUBMIT_DIR, filename);
                postReport(file);
                logger.info("Wrote report to " + file);
            } else {
                logger.info("Not writing report: sentinel " + REPORT_SENTINEL + " doesn't exist");
            }
        } catch ( Exception e ) {
            // we catch everything, and no matter what eat the error
            logger.warn("Received error while posting report.  GATK continuing on but no run report has been generated because: " + e.getMessage());
        }
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

        public ExceptionToXML(Throwable e) {
            message = e.getMessage();
            for (StackTraceElement element : e.getStackTrace()) {
                stackTrace.add(element.toString());
            }

            if ( e.getCause() != null ) {
                cause = new ExceptionToXML(e.getCause());
            }
        }
    }
}
