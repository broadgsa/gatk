package org.broadinstitute.sting.gatk.traversals;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.filters.CountingFilteringIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.MathUtils;

import java.util.Map;
import java.util.List;
import java.util.Iterator;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMRecord;

public abstract class TraversalEngine<M,T,WalkerType extends Walker<M,T>,ProviderType extends ShardDataProvider> {
    // Time in milliseconds since we initialized this engine
    private long startTime = -1;
    private long lastProgressPrintTime = -1;                // When was the last time we printed our progress?

    // How long can we go without printing some progress info?
    private final long MAX_PROGRESS_PRINT_TIME = 30 * 1000;        // 10 seconds in millisecs
    private final long N_RECORDS_TO_PRINT = 1000000;

    // Maximum number of reads to process before finishing
    protected long maximumIterations = -1;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(TraversalEngine.class);

    /**
     * set the max number of iterations
     * @param maximumIterations the number of iterations
     */
    public void setMaximumIterations(final int maximumIterations) {
        this.maximumIterations = maximumIterations;
    }

    /**
     * @param curTime (current runtime, in millisecs)
     *
     * @return true if the maximum interval (in millisecs) has passed since the last printing
     */
    private boolean maxElapsedIntervalForPrinting(final long curTime) {
        return (curTime - this.lastProgressPrintTime) > MAX_PROGRESS_PRINT_TIME;
    }

    /**
     * Forward request to printProgress
     *
     * @param type the TRAVERSAL_TYPE of the traversal
     * @param loc  the location
     */
    public void printProgress(final String type, GenomeLoc loc) {
        printProgress(false, type, loc);
    }

    /**
     * Utility routine that prints out process information (including timing) every N records or
     * every M seconds, for N and M set in global variables.
     *
     * @param mustPrint If true, will print out info, regardless of nRecords or time interval
     * @param type      String to print out describing our atomic traversal type ("read", "locus", etc)
     * @param loc       Current location
     */
    private void printProgress(boolean mustPrint, final String type, GenomeLoc loc) {
        final long nRecords = TraversalStatistics.nRecords;
        final long curTime = System.currentTimeMillis();
        final double elapsed = (curTime - startTime) / 1000.0;
        //System.out.printf("Cur = %d, last print = %d, elapsed=%.2f, nRecords=%d, met=%b%n", curTime, lastProgressPrintTime, elapsed, nRecords, maxElapsedIntervalForPrinting(curTime));

        if (mustPrint || nRecords == 1 || nRecords % N_RECORDS_TO_PRINT == 0 || maxElapsedIntervalForPrinting(curTime)) {
            this.lastProgressPrintTime = curTime;
            final double secsPer1MReads = (elapsed * 1000000.0) / nRecords;
            if (loc != null)
                logger.info(String.format("[PROGRESS] Traversed to %s, processing %,d %s in %.2f secs (%.2f secs per 1M %s)", loc, nRecords, type, elapsed, secsPer1MReads, type));
            else
                logger.info(String.format("[PROGRESS] Traversed %,d %s in %.2f secs (%.2f secs per 1M %s)", nRecords, type, elapsed, secsPer1MReads, type));
        }   
    }

    /**
     * A passthrough method so that subclasses can report which types of traversals they're using.
     *
     * @param sum Result of the computation.
     */
    public abstract void printOnTraversalDone(T sum);

    /**
     * Called after a traversal to print out information about the traversal process
     *
     * @param type describing this type of traversal
     * @param sum  The reduce result of the traversal
     */
    protected void printOnTraversalDone(final String type, T sum) {
        printProgress(true, type, null);
        final long curTime = System.currentTimeMillis();
        final double elapsed = (curTime - startTime) / 1000.0;

        // count up the number of skipped reads by summing over all filters
        long nSkippedReads = 0L;
        for ( long counts : TraversalStatistics.counter.values() )
            nSkippedReads += counts;

        logger.info(String.format("Total runtime %.2f secs, %.2f min, %.2f hours%n", elapsed, elapsed / 60, elapsed / 3600));
        logger.info(String.format("%d reads were filtered out during traversal out of %d total (%.2f%%)",
                nSkippedReads,
                TraversalStatistics.nReads,
                100.0 * MathUtils.ratio(nSkippedReads, TraversalStatistics.nReads)));
        for ( Map.Entry<Class, Long> filterCounts : TraversalStatistics.counter.entrySet() ) {
            long count = filterCounts.getValue();
            logger.info(String.format("  -> %d reads (%.2f%% of total) failing %s",
                    count, 100.0 * MathUtils.ratio(count, TraversalStatistics.nReads), Utils.getClassName(filterCounts.getKey())));
        }
    }

    /** Initialize the traversal engine.  After this point traversals can be run over the data */
    public void initialize() {
        lastProgressPrintTime = startTime = System.currentTimeMillis();
    }

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

    public static Iterator<SAMRecord> addMandatoryFilteringIterators(Iterator<SAMRecord> iter, List<SamRecordFilter> filters ) {
        for( SamRecordFilter filter : filters ) {
            //logger.debug("Adding filter " + filter.getClass());
            iter = new CountingFilteringIterator(iter,filter);
        }

        return new CountingFilteringIterator(iter); // special case to count all reads
    }


}
