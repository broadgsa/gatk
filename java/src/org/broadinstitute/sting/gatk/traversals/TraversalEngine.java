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

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.ReadMetrics;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

public abstract class TraversalEngine<M,T,WalkerType extends Walker<M,T>,ProviderType extends ShardDataProvider> {
    // Time in milliseconds since we initialized this engine
    private static final int HISTORY_WINDOW_SIZE = 50;

    private static class ProcessingHistory {
        double elapsedSeconds;
        long unitsProcessed;
        long bpProcessed;
        GenomeLoc loc;

        public ProcessingHistory(double elapsedSeconds, GenomeLoc loc, long unitsProcessed, long bpProcessed) {
            this.elapsedSeconds = elapsedSeconds;
            this.loc = loc;
            this.unitsProcessed = unitsProcessed;
            this.bpProcessed = bpProcessed;
        }

    }

    /**
     * Simple utility class that makes it convenient to print unit adjusted times
     */
    private static class MyTime {
        double t;           // in Seconds
        int precision;      // for format

        public MyTime(double t, int precision) {
            this.t = t;
            this.precision = precision;
        }

        public MyTime(double t) {
            this(t, 1);
        }

        /**
         * Instead of 10000 s, returns 2.8 hours
         * @return
         */
        public String toString() {
            double unitTime = t;
            String unit = "s";

            if ( t > 120 ) {
                unitTime = t / 60; // minutes
                unit = "m";

                if ( unitTime > 120 ) {
                    unitTime /= 60; // hours
                    unit = "h";

                    if ( unitTime > 100 ) {
                        unitTime /= 24; // days
                        unit = "d";

                        if ( unitTime > 20 ) {
                            unitTime /= 7; // days
                            unit = "w";
                        }
                    }
                }
            }

            return String.format("%6."+precision+"f %s", unitTime, unit);
        }
    }

    /** lock object to sure updates to history are consistent across threads */
    private static final Object lock = new Object();
    LinkedList<ProcessingHistory> history = new LinkedList<ProcessingHistory>();

    /** We use the SimpleTimer to time our run */
    private SimpleTimer timer = new SimpleTimer("Traversal");

    // How long can we go without printing some progress info?
    private long lastProgressPrintTime = -1;                       // When was the last time we printed progress log?
    private final long PROGRESS_PRINT_FREQUENCY = 10 * 1000;        // in seconds

    // for performance log
    private static final boolean PERFORMANCE_LOG_ENABLED = true;
    private PrintStream performanceLog = null;
    private long lastPerformanceLogPrintTime = -1;                   // When was the last time we printed to the performance log?
    private final long PERFORMANCE_LOG_PRINT_FREQUENCY = PROGRESS_PRINT_FREQUENCY;  // in seconds

    /** Size, in bp, of the area we are processing.  Updated once in the system in initial for performance reasons */
    long targetSize = -1;
    GenomeLocSortedSet targetIntervals = null;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(TraversalEngine.class);

    protected GenomeAnalysisEngine engine;

    // ----------------------------------------------------------------------------------------------------
    //
    // ABSTRACT METHODS
    //
    // ----------------------------------------------------------------------------------------------------
    /**
     * Gets the named traversal type associated with the given traversal.
     * @return A user-friendly name for the given traversal type.
     */
    protected abstract String getTraversalType();

    /**
     * this method must be implemented by all traversal engines
     *
     * @param walker       the walker to run with
     * @param dataProvider the data provider that generates data given the shard
     * @param sum          the accumulator
     *
     * @return an object of the reduce type
     */
    public abstract T traverse(WalkerType walker,
                               ProviderType dataProvider,
                               T sum);

    // ----------------------------------------------------------------------------------------------------
    //
    // Common timing routines
    //
    // ----------------------------------------------------------------------------------------------------
    /**
     * Initialize the traversal engine.  After this point traversals can be run over the data
     * @param engine GenomeAnalysisEngine for this traversal
     */
    public void initialize(GenomeAnalysisEngine engine) {
        if ( engine == null )
            throw new ReviewedStingException("BUG: GenomeAnalysisEngine cannot be null!");

        this.engine = engine;

        if ( PERFORMANCE_LOG_ENABLED && engine.getArguments() != null && engine.getArguments().performanceLog != null ) {
            try {
                performanceLog = new PrintStream(new FileOutputStream(engine.getArguments().performanceLog));
                List<String> pLogHeader = Arrays.asList("elapsed.time", "units.processed", "processing.speed", "bp.processed", "bp.speed", "genome.fraction.complete", "est.total.runtime", "est.time.remaining");
                performanceLog.println(Utils.join("\t", pLogHeader));
            } catch (FileNotFoundException e) {
                throw new UserException.CouldNotCreateOutputFile(engine.getArguments().performanceLog, e);
            }
        }

        // if we don't have any intervals defined, create intervals from the reference itself
        if ( this.engine.getIntervals() == null )
            targetIntervals = GenomeLocSortedSet.createSetFromSequenceDictionary(engine.getReferenceDataSource().getReference().getSequenceDictionary());
        else
            targetIntervals = this.engine.getIntervals();
        targetSize = targetIntervals.coveredSize();
    }

    /**
     * Should be called to indicate that we're going to process records and the timer should start ticking
     */
    public void startTimers() {
        timer.start();
        lastProgressPrintTime = timer.currentTime();
    }

    /**
     * @param curTime (current runtime, in millisecs)
     * @param lastPrintTime the last time we printed, in machine milliseconds
     * @param printFreq maximum permitted difference between last print and current times
     *
     * @return true if the maximum interval (in millisecs) has passed since the last printing
     */
    private boolean maxElapsedIntervalForPrinting(final long curTime, long lastPrintTime, long printFreq) {
        return (curTime - lastPrintTime) > printFreq;
    }

    /**
     * Forward request to printProgress
     *
     * @param shard the given shard currently being processed.
     * @param loc  the location
     */
    public void printProgress(Shard shard, GenomeLoc loc) {
        // A bypass is inserted here for unit testing.
        // TODO: print metrics outside of the traversal engine to more easily handle cumulative stats.
        ReadMetrics cumulativeMetrics = engine.getCumulativeMetrics() != null ? engine.getCumulativeMetrics().clone() : new ReadMetrics();
        cumulativeMetrics.incrementMetrics(shard.getReadMetrics());
        printProgress(loc, cumulativeMetrics, false);
    }

    /**
     * Utility routine that prints out process information (including timing) every N records or
     * every M seconds, for N and M set in global variables.
     *
     * @param loc       Current location, can be null if you are at the end of the traversal
     * @param metrics   Metrics of reads filtered in/out.
     * @param mustPrint If true, will print out info, regardless of nRecords or time interval
     */
    private void printProgress(GenomeLoc loc, ReadMetrics metrics, boolean mustPrint) {
        final long nRecords = metrics.getNumIterations();

        if ( nRecords == 1 ) {
            logger.info("[INITIALIZATION COMPLETE; TRAVERSAL STARTING]");
            logger.info(String.format("%15s processed.%s  runtime per.1M.%s completed total.runtime remaining",
                    "Location", getTraversalType(), getTraversalType()));

        }
        else {
            final long curTime = timer.currentTime();
            boolean printProgress = mustPrint || maxElapsedIntervalForPrinting(curTime, lastProgressPrintTime, PROGRESS_PRINT_FREQUENCY);
            boolean printLog = performanceLog != null && maxElapsedIntervalForPrinting(curTime, lastPerformanceLogPrintTime, PERFORMANCE_LOG_PRINT_FREQUENCY);

            if ( printProgress || printLog ) {
                lastProgressPrintTime = curTime;
                ProcessingHistory last = updateHistory(loc, metrics);

                final MyTime elapsed = new MyTime(last.elapsedSeconds);
                final MyTime bpRate = new MyTime(secondsPerMillionBP(last));
                final MyTime unitRate = new MyTime(secondsPerMillionElements(last));
                final double fractionGenomeTargetCompleted = calculateFractionGenomeTargetCompleted(last);
                final MyTime estTotalRuntime = new MyTime(elapsed.t / fractionGenomeTargetCompleted);
                final MyTime timeToCompletion = new MyTime(estTotalRuntime.t - elapsed.t);

                if ( printProgress ) {
//                    String common = String.format("%4.1e %s in %s, %s per 1M %s, %4.1f%% complete, est. runtime %s, %s remaining",
//                            nRecords*1.0, getTraversalType(), elapsed, unitRate,
//                            getTraversalType(), 100*fractionGenomeTargetCompleted,
//                            estTotalRuntime, timeToCompletion);
//
//                    if (loc != null)
//                        logger.info(String.format("%20s: processing %s", loc, common));
//                    else
//                        logger.info(String.format("Processing %s", common));
                    logger.info(String.format("%15s        %5.2e %s     %s     %4.1f%%      %s  %s",
                            loc == null ? "done" : loc, nRecords*1.0, elapsed, unitRate,
                            100*fractionGenomeTargetCompleted, estTotalRuntime, timeToCompletion));

                }

                if ( printLog ) {
                    lastPerformanceLogPrintTime = curTime;
                    performanceLog.printf("%.2f\t%d\t%.2e\t%d\t%.2e\t%.2e\t%.2f\t%.2f%n",
                            elapsed.t, nRecords, unitRate.t, last.bpProcessed, bpRate.t,
                            fractionGenomeTargetCompleted, estTotalRuntime.t, timeToCompletion.t);
                }
            }
        }
    }

    /**
     * Keeps track of the last HISTORY_WINDOW_SIZE data points for the progress meter.  Currently the
     * history isn't used in any way, but in the future it'll become valuable for more accurate estimates
     * for when a process will complete.
     *
     * @param loc our current position.  If null, assumes we are done traversing
     * @param metrics information about what's been processed already
     * @return
     */
    private final ProcessingHistory updateHistory(GenomeLoc loc, ReadMetrics metrics) {
        synchronized (lock) {
            if ( history.size() > HISTORY_WINDOW_SIZE )
                history.pop();

            long nRecords = metrics.getNumIterations();
            long bpProcessed = loc == null ? targetSize : targetIntervals.sizeBeforeLoc(loc); // null -> end of processing
            history.add(new ProcessingHistory(timer.getElapsedTime(), loc, nRecords, bpProcessed));

            return history.getLast();
        }
    }

    /** How long in seconds to process 1M traversal units? */
    private final double secondsPerMillionElements(ProcessingHistory last) {
        return (last.elapsedSeconds * 1000000.0) / Math.max(last.unitsProcessed, 1);
    }

    /** How long in seconds to process 1M bp on the genome? */
    private final double secondsPerMillionBP(ProcessingHistory last) {
        return (last.elapsedSeconds * 1000000.0) / Math.max(last.bpProcessed, 1);
    }

    /** What fractoin of the target intervals have we covered? */
    private final double calculateFractionGenomeTargetCompleted(ProcessingHistory last) {
        return (1.0*last.bpProcessed) / targetSize;
    }

    /**
     * Called after a traversal to print out information about the traversal process
     */
    public void printOnTraversalDone(ReadMetrics cumulativeMetrics) {
        printProgress(null, cumulativeMetrics, true);

        final double elapsed = timer.getElapsedTime();

        // count up the number of skipped reads by summing over all filters
        long nSkippedReads = 0L;
        for ( Map.Entry<Class, Long> countsByFilter: cumulativeMetrics.getCountsByFilter().entrySet())
            nSkippedReads += countsByFilter.getValue();

        logger.info(String.format("Total runtime %.2f secs, %.2f min, %.2f hours", elapsed, elapsed / 60, elapsed / 3600));
        if ( cumulativeMetrics.getNumReadsSeen() > 0 )
            logger.info(String.format("%d reads were filtered out during traversal out of %d total (%.2f%%)",
                    nSkippedReads,
                    cumulativeMetrics.getNumReadsSeen(),
                    100.0 * MathUtils.ratio(nSkippedReads,cumulativeMetrics.getNumReadsSeen())));
        for ( Map.Entry<Class, Long> filterCounts : cumulativeMetrics.getCountsByFilter().entrySet() ) {
            long count = filterCounts.getValue();
            logger.info(String.format("  -> %d reads (%.2f%% of total) failing %s",
                    count, 100.0 * MathUtils.ratio(count,cumulativeMetrics.getNumReadsSeen()), Utils.getClassName(filterCounts.getKey())));
        }

        if ( performanceLog != null ) performanceLog.close();
    }
}
