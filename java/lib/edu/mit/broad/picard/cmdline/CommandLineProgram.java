/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.cmdline;

import edu.mit.broad.picard.util.Log;
import edu.mit.broad.picard.util.StringUtil;
import edu.mit.broad.picard.metrics.Header;
import edu.mit.broad.picard.metrics.StringHeader;
import edu.mit.broad.picard.metrics.MetricsFile;
import edu.mit.broad.picard.metrics.MetricBase;

import java.io.File;
import java.util.Date;
import java.util.List;
import java.util.ArrayList;

/**
 * Abstract class to facilitate writing command-line programs.
 *
 * To use:
 *
 * 1. Extend this class with a concrete class that has data members annotated with @Option, @PositionalArguments
 * and/or @Usage annotations.
 *
 * 2. If there is any custom command-line validation, override customCommandLineValidation().  When this method is
 * called, the command line has been parsed and set into the data members of the concrete class.
 *
 * 3. Implement a method doWork().  This is called after successful comand-line processing.  The value it returns is
 * the exit status of the program.  It is assumed that the concrete class emits any appropriate error message before
 * returning non-zero.  doWork() may throw unchecked exceptions, which are caught and reported appropriately.
 *
 * 4. Implement the following static method in the concrete class:
 *
 *     public static void main(String[] argv) {
        System.exit(new MyConcreteClass().instanceMain(argv));
    }


 */
public abstract class CommandLineProgram {

    @Option
    public File TMP_DIR = new File(System.getProperty("java.io.tmpdir"), System.getProperty("user.name"));

    @Option(doc = "Control verbosity of logging")
    public Log.LogLevel VERBOSITY = Log.LogLevel.INFO;

    @Option(doc = "Whether to suppress job-summary info on System.out")
    public Boolean QUIET = false;

    private final String standardUsagePreamble = CommandLineParser.getStandardUsagePreamble(getClass());

    /**
     * Initialized in parseArgs.  Subclasses may want to access this to do
     * their own validation, and then print usage using clp.
     */
    protected CommandLineParser clp;

    private final List<Header> defaultHeaders = new ArrayList<Header>();

    /**
     * Do the work after command line has been parsed.
     * RuntimeException may be thrown by this method, and are reported appropriately.
     * @return program exit status.
     */
    protected abstract int doWork();

    public int instanceMain(final String[] argv) {
        // Build the default headers
        final Date startDate = new Date();
        final String cmdline = getClass().getName() + " " + StringUtil.join(" ", argv);
        this.defaultHeaders.add(new StringHeader(cmdline));
        this.defaultHeaders.add(new StringHeader("Started on: " + startDate));

        if (!parseArgs(argv)) {
            return 1;
        }

        Log.setGlobalLogLevel(VERBOSITY);

        if (!TMP_DIR.exists()) {
            // Intentially not checking the return value, because it may be that the program does not
            // need a tmp_dir.  If this fails, the problem will be discovered downstream.
            TMP_DIR.mkdir();
        }
        System.setProperty("java.io.tmpdir", TMP_DIR.getAbsolutePath());
        if (!QUIET) {
            System.out.println("[" + new Date() + "] " + cmdline);
        }
        final int ret = doWork();
        if (!QUIET) {
            System.out.println("[" + new Date() + "] " + getClass().getName() + " done.");
            System.out.println("Runtime.totalMemory()=" + Runtime.getRuntime().totalMemory());
        }
        return ret;
    }

    /**
     * Put any custom command-line validation in an override of this method.
     * clp is initialized at this point and can be used to print usage and access argv.
     * Any options set by command-line parser can be validated.
     * @return true if command line is valid.
     */
    protected boolean customCommandLineValidation() {
        return true;
    }

    /**
     *
     * @return true if command line is valid
     */
    protected boolean parseArgs(final String[] argv) {
        clp = new CommandLineParser(this);
        final boolean ret = clp.parseOptions(System.err, argv);
        if (!ret) {
            return false;
        }
        return customCommandLineValidation();
    }

    /** Gets a MetricsFile with default headers already written into it. */
    protected <A extends MetricBase,B extends Comparable> MetricsFile<A,B> getMetricsFile() {
        final MetricsFile<A,B> file = new MetricsFile<A,B>();
        for (final Header h : this.defaultHeaders) {
            file.addHeader(h);
        }

        return file;
    }

    public String getStandardUsagePreamble() {
        return standardUsagePreamble;
    }
}
