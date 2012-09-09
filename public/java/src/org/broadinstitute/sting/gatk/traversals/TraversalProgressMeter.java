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

package org.broadinstitute.sting.gatk.traversals;

import com.google.java.contract.Ensures;
import com.google.java.contract.Invariant;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.ReadMetrics;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * A progress meter that prints a few key metrics to a logger and optionally to a file
 *
 * Metrics include:
 *      -- Number of processed X (X = traversal units)
 *      -- Runtime per.1M X
 *      -- Percent of regions to be processed completed
 *      -- The estimated total runtime based on previous performance
 *      -- The estimated time remaining for the entire process
 *
 * The optional file log an expanded set of metrics in tabular format
 * suitable for subsequent analysis in R.
 *
 * This class is -- and MUST BE -- thread-safe for use in the GATK.  Multiple independent
 * threads executing traversals will be calling printProgress() simultaneously and this
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
public class TraversalProgressMeter {
    protected static final Logger logger = Logger.getLogger(TraversalProgressMeter.class);

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
     * Used to get the total number of records we've processed so far.
     */
    final ReadMetrics cumulativeMetrics;

    /**
     * A string describing the type of this traversal, so we can say things like
     * "we are running at X traversalType per second"
     */
    private final String traversalType;

    /**
     * A potentially null file where we print a supplementary, R readable performance log
     * file.
     */
    private final PrintStream performanceLog;

    /** We use the SimpleTimer to time our run */
    private final SimpleTimer timer = new SimpleTimer("Traversal");

    /**
     * Create a new TraversalProgressMeter
     *
     * @param cumulativeMetrics the object where the shared traversal counts are being updated
     * @param performanceLogFile an optional performance log file where a table of performance logs will be written
     * @param traversalUnits the name of this traversal type, suitable for saying X seconds per traversalUnits
     * @param processingIntervals the intervals being processed
     */
    public TraversalProgressMeter(final ReadMetrics cumulativeMetrics,
                                  final File performanceLogFile,
                                  final String traversalUnits,
                                  final GenomeLocSortedSet processingIntervals) {
        if ( cumulativeMetrics == null ) throw new IllegalArgumentException("cumulativeMetrics cannot be null!");
        if ( traversalUnits == null ) throw new IllegalArgumentException("traversalUnits cannot be null");
        if ( processingIntervals == null ) throw new IllegalArgumentException("Target intervals cannot be null");

        this.cumulativeMetrics = cumulativeMetrics;
        this.traversalType = traversalUnits;
        this.regionsBeingProcessed = processingIntervals;

        // setup the performance logger output, if requested by the GATK engine
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
        start();
    }

    /**
     * Forward request to printProgress
     *
     * Assumes that one cycle has been completed
     *
     * @param loc  the location
     */
    public void printProgress(final GenomeLoc loc) {
        // A bypass is inserted here for unit testing.
        printProgress(loc, false);
    }

    private synchronized void start() {
        timer.start();
        lastProgressPrintTime = timer.currentTime();

        logger.info("[INITIALIZATION COMPLETE; TRAVERSAL STARTING]");
        logger.info(String.format("%15s processed.%s  runtime per.1M.%s completed total.runtime remaining",
                "Location", traversalType, traversalType));
    }

    /**
     * Utility routine that prints out process information (including timing) every N records or
     * every M seconds, for N and M set in global variables.
     *
     * Synchronized to ensure that even with multiple threads calling printProgress we still
     * get one clean stream of meter logs.
     *
     * @param loc       Current location, can be null if you are at the end of the traversal
     * @param mustPrint If true, will print out info, regardless of time interval
     */
    private synchronized void printProgress(final GenomeLoc loc, boolean mustPrint) {
        final long curTime = timer.currentTime();
        final boolean printProgress = mustPrint || maxElapsedIntervalForPrinting(curTime, lastProgressPrintTime, progressPrintFrequency);
        final boolean printLog = performanceLog != null && maxElapsedIntervalForPrinting(curTime, lastPerformanceLogPrintTime, PERFORMANCE_LOG_PRINT_FREQUENCY);

        if ( printProgress || printLog ) {
            final ProgressData progressData = takeProgressSnapshot(loc, cumulativeMetrics);

            final AutoFormattingTime elapsed = new AutoFormattingTime(progressData.elapsedSeconds);
            final AutoFormattingTime bpRate = new AutoFormattingTime(progressData.secondsPerMillionBP());
            final AutoFormattingTime unitRate = new AutoFormattingTime(progressData.secondsPerMillionElements());
            final double fractionGenomeTargetCompleted = progressData.calculateFractionGenomeTargetCompleted(targetSizeInBP);
            final AutoFormattingTime estTotalRuntime = new AutoFormattingTime(elapsed.getTimeInSeconds() / fractionGenomeTargetCompleted);
            final AutoFormattingTime timeToCompletion = new AutoFormattingTime(estTotalRuntime.getTimeInSeconds() - elapsed.getTimeInSeconds());

            if ( printProgress ) {
                lastProgressPrintTime = curTime;
                updateLoggerPrintFrequency(estTotalRuntime.getTimeInSeconds());

                // a pretty name for our position
                final String posName = loc == null
                        ? (mustPrint ? "done" : "unmapped reads")
                        : String.format("%s:%d", loc.getContig(), loc.getStart());

                logger.info(String.format("%15s        %5.2e %s     %s    %5.1f%%      %s  %s",
                        posName, progressData.unitsProcessed*1.0, elapsed, unitRate,
                        100*fractionGenomeTargetCompleted, estTotalRuntime, timeToCompletion));

            }

            if ( printLog ) {
                lastPerformanceLogPrintTime = curTime;
                performanceLog.printf("%.2f\t%d\t%.2e\t%d\t%.2e\t%.2e\t%.2f\t%.2f%n",
                        elapsed.getTimeInSeconds(), progressData.unitsProcessed, unitRate.getTimeInSeconds(),
                        progressData.bpProcessed, bpRate.getTimeInSeconds(),
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
     * @param metrics information about what's been processed already
     * @return
     */
    private ProgressData takeProgressSnapshot(final GenomeLoc loc, final ReadMetrics metrics) {
        final long nRecords = metrics.getNumIterations();
        // null -> end of processing
        final long bpProcessed = loc == null ? targetSizeInBP : regionsBeingProcessed.sizeBeforeLoc(loc);
        return new ProgressData(timer.getElapsedTime(), nRecords, bpProcessed);
    }

    /**
     * Called after a traversal to print out information about the traversal process
     */
    public void printOnDone() {
        printProgress(null, true);

        final double elapsed = timer == null ? 0 : timer.getElapsedTime();

        // count up the number of skipped reads by summing over all filters
        long nSkippedReads = 0L;
        for ( final long countsByFilter : cumulativeMetrics.getCountsByFilter().values())
            nSkippedReads += countsByFilter;

        logger.info(String.format("Total runtime %.2f secs, %.2f min, %.2f hours", elapsed, elapsed / 60, elapsed / 3600));

        // TODO -- move into MicroScheduler
        if ( cumulativeMetrics.getNumReadsSeen() > 0 )
            logger.info(String.format("%d reads were filtered out during traversal out of %d total (%.2f%%)",
                    nSkippedReads,
                    cumulativeMetrics.getNumReadsSeen(),
                    100.0 * MathUtils.ratio(nSkippedReads,cumulativeMetrics.getNumReadsSeen())));
        for ( Map.Entry<String, Long> filterCounts : cumulativeMetrics.getCountsByFilter().entrySet() ) {
            long count = filterCounts.getValue();
            logger.info(String.format("  -> %d reads (%.2f%% of total) failing %s",
                    count, 100.0 * MathUtils.ratio(count,cumulativeMetrics.getNumReadsSeen()), filterCounts.getKey()));
        }

        if ( performanceLog != null ) performanceLog.close();
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

    /**
     * a snapshot of our performance, suitable for storage and later analysis
     */
    private static class ProgressData {
        final double elapsedSeconds;
        final long unitsProcessed;
        final long bpProcessed;

        @Requires({"unitsProcessed >= 0", "bpProcessed >= 0", "elapsedSeconds >= 0"})
        public ProgressData(double elapsedSeconds, long unitsProcessed, long bpProcessed) {
            this.elapsedSeconds = elapsedSeconds;
            this.unitsProcessed = unitsProcessed;
            this.bpProcessed = bpProcessed;
        }

        /** How long in seconds to process 1M traversal units? */
        @Ensures("result >= 0.0")
        private double secondsPerMillionElements() {
            return (elapsedSeconds * 1000000.0) / Math.max(unitsProcessed, 1);
        }

        /** How long in seconds to process 1M bp on the genome? */
        @Ensures("result >= 0.0")
        private double secondsPerMillionBP() {
            return (elapsedSeconds * 1000000.0) / Math.max(bpProcessed, 1);
        }

        /** What fraction of the target intervals have we covered? */
        @Requires("targetSize >= 0")
        @Ensures({"result >= 0.0", "result <= 1.0"})
        private double calculateFractionGenomeTargetCompleted(final long targetSize) {
            return (1.0*bpProcessed) / Math.max(targetSize, 1);
        }
    }
}
