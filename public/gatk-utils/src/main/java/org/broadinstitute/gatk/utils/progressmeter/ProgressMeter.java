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

package org.broadinstitute.gatk.utils.progressmeter;

import com.google.java.contract.Ensures;
import com.google.java.contract.Invariant;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.*;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

/**
 * A meter measuring progress on a calculation through a set of genomic regions that can
 * print a few key metrics to a logger and optionally to a file
 *
 * The key information for assessing progress is a set of genome locs describing the total
 * set of regions we will process.  Whenever (at reasonable intervals) the processing unit
 * can called notifyOfProgress and this logger may, depending on the metering delay, print
 * a log message with the following metrics:
 *
 *      -- Number of processed X (X = processing units)
 *      -- Runtime per.1M X
 *      -- Percent of regions to be processed completed
 *      -- The estimated total runtime based on previous performance
 *      -- The estimated time remaining for the entire process
 *
 * The optional file log an expanded set of metrics in tabular format
 * suitable for subsequent analysis in R.
 *
 * This class is -- and MUST BE -- thread-safe for use in the GATK.  Multiple independent
 * threads executing processors will be calling notifyOfProgress() simultaneously and this
 * class does (and MUST) properly sort out the timings of logs without interlacing outputs
 * because of these threads.
 *
 * Consequently, the fundamental model for when to print the logs is time based.  We basically
 * print a meter message every X seconds, minutes, hours, whatever is appropriate based on the
 * estimated remaining runtime.
 *
 * @author depristo
 * @since 2010 maybe, but written in 09/12 for clarity
 */
@Invariant({
        "targetSizeInBP >= 0",
        "progressPrintFrequency > 0"
})
public class ProgressMeter {
    protected static final Logger logger = Logger.getLogger(ProgressMeter.class);

    // --------------------------------------------------------------------------------
    // static constants controlling overall system behavior
    // --------------------------------------------------------------------------------

    /**
     * Min. milliseconds after we start up the meter before we will print our first meter message
     */
    private final static long MIN_ELAPSED_TIME_BEFORE_FIRST_PROGRESS = 30 * 1000;

    /**
     * How often should we print performance logging information, when we are sending this
     * information to a file?  Not dynamically updated as the logger meter is.
     */
    private final static long PERFORMANCE_LOG_PRINT_FREQUENCY = 10 * 1000;

    private final static double TWO_HOURS_IN_SECONDS    =  2.0 * 60.0 * 60.0;
    private final static double TWELVE_HOURS_IN_SECONDS = 12.0 * 60.0 * 60.0;

    // --------------------------------------------------------------------------------
    // Variables we updating during running
    // --------------------------------------------------------------------------------

    /**
     * When was the last time we printed progress log?  In milleseconds
     */
    private long lastProgressPrintTime = -1;

    /**
     * How frequently should we be printing our meter messages?  Dynamically updated
     * depending on how long we think the run has left.
     */
    private long progressPrintFrequency = 10 * 1000; // default value

    /**
     * When was the last time we printed to the performance log?  In millseconds
     */
    private long lastPerformanceLogPrintTime = -1;

    // --------------------------------------------------------------------------------
    // final variables fixed at object creation time
    // --------------------------------------------------------------------------------

    /**
     * The set of genome locs describing the total region we are processing with
     * this GATK run.  Used to determine how close we are to completing the run
     */
    private final GenomeLocSortedSet regionsBeingProcessed;

    /**
     * Size, in bp, of the area we are processing, derived from regionsBeingProcessed.
     * Updated once in the system in initial for performance reasons
     */
    private final long targetSizeInBP;

    /**
     * A string describing the type of units being processes, so we can say things like
     * "we are running at X processingUnitName per second"
     */
    private final String processingUnitName;

    /**
     * The space allocated to #processingUnitName in the output
     */
    private final int processingUnitWidth;

    /**
     * The format string used for progress lines
     */
    private final String progressFormatString;

    /**
     * A potentially null file where we print a supplementary, R readable performance log
     * file.
     */
    private final PrintStream performanceLog;

    /** We use the SimpleTimer to time our run */
    private final SimpleTimer timer = new SimpleTimer();

    private GenomeLoc maxGenomeLoc = null;
    private Position position = new Position(PositionStatus.STARTING);
    private long nTotalRecordsProcessed = 0;

    /**
     * The elapsed time in nanosecond, updated by the daemon thread, so that
     * we don't pay any system call overhead to determine the the elapsed time.
     */
    private long elapsedTimeInNanosecondUpdatedByDaemon = 0;

    final ProgressMeterDaemon progressMeterDaemon;

    /**
     * Create a new ProgressMeter
     *
     * Note that progress meter isn't started until the client calls start()
     *
     * @param performanceLogFile an optional performance log file where a table of performance logs will be written
     * @param processingUnitName the name of the unit type being processed, suitable for saying X seconds per processingUnitName
     * @param processingIntervals the intervals being processed
     */
    public ProgressMeter(final File performanceLogFile,
                         final String processingUnitName,
                         final GenomeLocSortedSet processingIntervals) {
        this(performanceLogFile, processingUnitName, processingIntervals, ProgressMeterDaemon.DEFAULT_POLL_FREQUENCY_MILLISECONDS);
    }

    protected ProgressMeter(final File performanceLogFile,
                            final String processingUnitName,
                            final GenomeLocSortedSet processingIntervals,
                            final long pollingFrequency) {
        if ( processingUnitName == null ) throw new IllegalArgumentException("processingUnitName cannot be null");
        if ( processingIntervals == null ) throw new IllegalArgumentException("Target intervals cannot be null");

        this.processingUnitName = processingUnitName;
        this.regionsBeingProcessed = processingIntervals;
        this.processingUnitWidth = Math.max(processingUnitName.length(), "processed".length());
        this.progressFormatString = String.format("%%15s   %%%1$ds   %%7s   %%%1$ds      %%5.1f%%%%   %%7s   %%9s", processingUnitWidth);

        // setup the performance logger output, if requested
        if ( performanceLogFile != null ) {
            try {
                this.performanceLog = new PrintStream(new FileOutputStream(performanceLogFile));
                final List<String> pLogHeader = Arrays.asList("elapsed.time", "units.processed", "processing.speed",
                        "bp.processed", "bp.speed", "genome.fraction.complete", "est.total.runtime", "est.time.remaining");
                performanceLog.println(Utils.join("\t", pLogHeader));
            } catch (FileNotFoundException e) {
                throw new UserException.CouldNotCreateOutputFile(performanceLogFile, e);
            }
        } else {
            performanceLog = null;
        }

        // cached for performance reasons
        targetSizeInBP = processingIntervals.coveredSize();

        // start up the timer
        progressMeterDaemon = new ProgressMeterDaemon(this, pollingFrequency);
    }

    public ProgressMeterDaemon getProgressMeterDaemon() {
        return progressMeterDaemon;
    }

    /**
     * Start up the progress meter, printing initialization message and starting up the
     * daemon thread for periodic printing.
     */
    @Requires("progressMeterDaemon != null")
    public synchronized void start() {
        timer.start();
        lastProgressPrintTime = timer.currentTime();
        final String formatString = String.format("%%15s | %%%1$ds | %%7s | %%%1$ds | %%9s | %%7s | %%9s", processingUnitWidth);

        logger.info("[INITIALIZATION COMPLETE; STARTING PROCESSING]");
        logger.info(String.format(formatString, "", "processed", "time", "per 1M", "", "total", "remaining"));
        logger.info(String.format(formatString, "Location", processingUnitName, "elapsed", processingUnitName,
                "completed", "runtime", "runtime"));

        progressMeterDaemon.start();
    }

    /**
     * @return the current runtime in nanoseconds
     */
    @Ensures("result >= 0")
    public long getRuntimeInNanoseconds() {
        return timer.getElapsedTimeNano();
    }

    /**
     * This function is just like getRuntimeInNanoseconds but it doesn't actually query the
     * system timer to determine the value, but rather uses a local variable in this meter
     * that is updated by the daemon thread.  This means that the result is ridiculously imprecise
     * for a nanosecond value (as it's only updated each pollingFrequency of the daemon) but
     * it is free for clients to access, which can be critical when one wants to do tests like:
     *
     * for some work unit:
     *   do unit if getRuntimeInNanosecondsUpdatedPeriodically < X
     *
     * and have this operation eventually timeout but don't want to pay the system call time to
     * ensure that the loop exits as soon as the elapsed time exceeds X
     *
     * @return the current runtime in nanoseconds
     */
    @Ensures("result >= 0")
    public long getRuntimeInNanosecondsUpdatedPeriodically() {
        return elapsedTimeInNanosecondUpdatedByDaemon;
    }

    /**
     * Update the period runtime variable to the current runtime in nanoseconds.  Should only
     * be called by the daemon thread
     */
    protected void updateElapsedTimeInNanoseconds() {
        elapsedTimeInNanosecondUpdatedByDaemon = getRuntimeInNanoseconds();
    }



    /**
     * Utility routine that prints out process information (including timing) every N records or
     * every M seconds, for N and M set in global variables.
     *
     * Synchronized to ensure that even with multiple threads calling notifyOfProgress we still
     * get one clean stream of meter logs.
     *
     * Note this thread doesn't actually print progress, unless must print is true, but just registers
     * the progress itself.  A separate printing daemon periodically polls the meter to print out
     * progress
     *
     * @param loc Current location, can be null if you are at the end of the processing unit.  Must
     *            have size == 1 (cannot be multiple bases in size).
     * @param nTotalRecordsProcessed the total number of records we've processed
     */
    public synchronized void notifyOfProgress(final GenomeLoc loc, final long nTotalRecordsProcessed) {
        if ( nTotalRecordsProcessed < 0 ) throw new IllegalArgumentException("nTotalRecordsProcessed must be >= 0");
        if ( loc.size() != 1 ) throw new IllegalArgumentException("GenomeLoc must have size == 1 but got " + loc);

        // weird comparison to ensure that loc == null (in unmapped reads) is keep before maxGenomeLoc == null (on startup)
        this.maxGenomeLoc = loc == null ? loc : (maxGenomeLoc == null ? loc : loc.max(maxGenomeLoc));
        this.nTotalRecordsProcessed = Math.max(this.nTotalRecordsProcessed, nTotalRecordsProcessed);

        // a pretty name for our position
        this.position = maxGenomeLoc == null ? new Position(PositionStatus.IN_UNMAPPED_READS) : new Position(maxGenomeLoc);
    }

    /**
     * Describes the status of this position marker, such as starting up, done, in the unmapped reads,
     * or somewhere on the genome
     */
    private enum PositionStatus {
        STARTING("Starting"),
        DONE("done"),
        IN_UNMAPPED_READS("unmapped reads"),
        ON_GENOME(null);

        public final String message;

        private PositionStatus(String message) {
            this.message = message;
        }
    }

    /**
     * A pair of position status and the genome loc, if necessary.  Used to get a
     * status update message as needed, without the computational cost of formatting
     * the genome loc string every time a progress notification happens (which is almost
     * always not printed)
     */
    private class Position {
        final PositionStatus type;
        final GenomeLoc maybeLoc;

        /**
         * Create a position object of any type != ON_GENOME
         * @param type
         */
        @Requires({"type != null", "type != PositionStatus.ON_GENOME"})
        private Position(PositionStatus type) {
            this.type = type;
            this.maybeLoc = null;
        }

        /**
         * Create a position object of type ON_GENOME at genomeloc loc
         * @param loc
         */
        @Requires("loc != null")
        private Position(GenomeLoc loc) {
            this.type = PositionStatus.ON_GENOME;
            this.maybeLoc = loc;
        }

        /**
         * @return a human-readable representation of this position
         */
        private String getMessage() {
            if ( type == PositionStatus.ON_GENOME )
                return maxGenomeLoc.getContig() + ":" + maxGenomeLoc.getStart();
            else
                return type.message;
        }
    }

    /**
     * Actually try to print out progress
     *
     * This function may print out if the progress print is due, but if not enough time has elapsed
     * since the last print we will not print out information.
     *
     * @param mustPrint if true, progress will be printed regardless of the last time we printed progress
     */
    protected synchronized void printProgress(final boolean mustPrint) {
        final long curTime = timer.currentTime();
        final boolean printProgress = mustPrint || maxElapsedIntervalForPrinting(curTime, lastProgressPrintTime, progressPrintFrequency);
        final boolean printLog = performanceLog != null && maxElapsedIntervalForPrinting(curTime, lastPerformanceLogPrintTime, PERFORMANCE_LOG_PRINT_FREQUENCY);

        if ( printProgress || printLog ) {
            final ProgressMeterData progressData = takeProgressSnapshot(maxGenomeLoc, nTotalRecordsProcessed);

            final AutoFormattingTime elapsed = new AutoFormattingTime(progressData.getElapsedSeconds(), 5, 1);
            final AutoFormattingTime bpRate = new AutoFormattingTime(progressData.secondsPerMillionBP());
            final AutoFormattingTime unitRate = new AutoFormattingTime(progressData.secondsPerMillionElements());
            final double fractionGenomeTargetCompleted = progressData.calculateFractionGenomeTargetCompleted(targetSizeInBP);
            final AutoFormattingTime estTotalRuntime = new AutoFormattingTime(elapsed.getTimeInSeconds() / fractionGenomeTargetCompleted, 5, 1);
            final AutoFormattingTime timeToCompletion = new AutoFormattingTime(estTotalRuntime.getTimeInSeconds() - elapsed.getTimeInSeconds());

            if ( printProgress ) {
                lastProgressPrintTime = curTime;
                updateLoggerPrintFrequency(estTotalRuntime.getTimeInSeconds());

                logger.info(String.format(progressFormatString,
                        position.getMessage(), progressData.getUnitsProcessed()*1.0, elapsed, unitRate,
                        100*fractionGenomeTargetCompleted, estTotalRuntime, timeToCompletion));

            }

            if ( printLog ) {
                lastPerformanceLogPrintTime = curTime;
                performanceLog.printf("%.2f\t%d\t%.2e\t%d\t%.2e\t%.2e\t%.2f\t%.2f%n",
                        elapsed.getTimeInSeconds(), progressData.getUnitsProcessed(), unitRate.getTimeInSeconds(),
                        progressData.getBpProcessed(), bpRate.getTimeInSeconds(),
                        fractionGenomeTargetCompleted, estTotalRuntime.getTimeInSeconds(),
                        timeToCompletion.getTimeInSeconds());
            }
        }
    }

    /**
     * Determine, based on remaining runtime, how often to print the meter
     *
     * @param totalRuntimeSeconds kinda obvious, no?
     */
    private void updateLoggerPrintFrequency(final double totalRuntimeSeconds) {
        // dynamically change the update rate so that short running jobs receive frequent updates while longer jobs receive fewer updates
        if ( totalRuntimeSeconds > TWELVE_HOURS_IN_SECONDS )
            progressPrintFrequency = 60 * 1000; // in milliseconds
        else if ( totalRuntimeSeconds > TWO_HOURS_IN_SECONDS )
            progressPrintFrequency = 30 * 1000; // in milliseconds
        else
            progressPrintFrequency = 10 * 1000; // in milliseconds
    }

    /**
     * Creates a new ProgressData object recording a snapshot of our progress at this instant
     *
     * @param loc our current position.  If null, assumes we are done traversing
     * @param nTotalRecordsProcessed the total number of records we've processed
     * @return
     */
    private ProgressMeterData takeProgressSnapshot(final GenomeLoc loc, final long nTotalRecordsProcessed) {
        // null -> end of processing
        final long bpProcessed = loc == null ? targetSizeInBP : regionsBeingProcessed.sizeBeforeLoc(loc);
        return new ProgressMeterData(timer.getElapsedTime(), nTotalRecordsProcessed, bpProcessed);
    }

    /**
     * Should be called when processing is done
     */
    public void notifyDone(final long nTotalRecordsProcessed) {
        // print out the progress meter
        this.nTotalRecordsProcessed = nTotalRecordsProcessed;
        this.position = new Position(PositionStatus.DONE);
        printProgress(true);

        logger.info(String.format("Total runtime %.2f secs, %.2f min, %.2f hours",
                timer.getElapsedTime(), timer.getElapsedTime() / 60, timer.getElapsedTime() / 3600));

        if ( performanceLog != null )
            performanceLog.close();

        // shutdown our daemon thread
        progressMeterDaemon.done();
    }

    /**
     * @param curTime (current runtime, in millisecs)
     * @param lastPrintTime the last time we printed, in machine milliseconds
     * @param printFreq maximum permitted difference between last print and current times
     *
     * @return true if the maximum interval (in millisecs) has passed since the last printing
     */
    private boolean maxElapsedIntervalForPrinting(final long curTime, long lastPrintTime, long printFreq) {
        final long elapsed = curTime - lastPrintTime;
        return elapsed > printFreq && elapsed > MIN_ELAPSED_TIME_BEFORE_FIRST_PROGRESS;
    }
}
