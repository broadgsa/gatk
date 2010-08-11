package org.broadinstitute.sting.gatk.traversals;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.ReadMetrics;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.MathUtils;

import java.util.Map;

public abstract class TraversalEngine<M,T,WalkerType extends Walker<M,T>,ProviderType extends ShardDataProvider> {
    // Time in milliseconds since we initialized this engine
    private long startTime = -1;
    private long lastProgressPrintTime = -1;                // When was the last time we printed our progress?

    // How long can we go without printing some progress info?
    private final long MAX_PROGRESS_PRINT_TIME = 30 * 1000;        // 10 seconds in millisecs
    private final long N_RECORDS_TO_PRINT = 1000000;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(TraversalEngine.class);

    /**
     * Gets the named traversal type associated with the given traversal.
     * @return A user-friendly name for the given traversal type.
     */
    protected abstract String getTraversalType();

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
     * @param shard the given shard currently being processed.
     * @param loc  the location
     */
    public void printProgress(Shard shard,GenomeLoc loc) {
        // A bypass is inserted here for unit testing.
        // TODO: print metrics outside of the traversal engine to more easily handle cumulative stats.
        ReadMetrics cumulativeMetrics = GenomeAnalysisEngine.instance != null ? GenomeAnalysisEngine.instance.getCumulativeMetrics().clone() : new ReadMetrics();
        cumulativeMetrics.incrementMetrics(shard.getReadMetrics());
        printProgress(loc, cumulativeMetrics, false);
    }

    /**
     * Utility routine that prints out process information (including timing) every N records or
     * every M seconds, for N and M set in global variables.
     *
     * @param loc       Current location
     * @param metrics   Metrics of reads filtered in/out.
     * @param mustPrint If true, will print out info, regardless of nRecords or time interval
     */
    private void printProgress(GenomeLoc loc, ReadMetrics metrics, boolean mustPrint) {
        final long nRecords = metrics.getNumIterations();
        final long curTime = System.currentTimeMillis();
        final double elapsed = (curTime - startTime) / 1000.0;
        //System.out.printf("Cur = %d, last print = %d, elapsed=%.2f, nRecords=%d, met=%b%n", curTime, lastProgressPrintTime, elapsed, nRecords, maxElapsedIntervalForPrinting(curTime));

        if (mustPrint || nRecords == 1 || nRecords % N_RECORDS_TO_PRINT == 0 || maxElapsedIntervalForPrinting(curTime)) {
            this.lastProgressPrintTime = curTime;
            final double secsPer1MReads = (elapsed * 1000000.0) / nRecords;
            if (loc != null)
                logger.info(String.format("[PROGRESS] Traversed to %s, processing %,d %s in %.2f secs (%.2f secs per 1M %s)", loc, nRecords, getTraversalType(), elapsed, secsPer1MReads, getTraversalType()));
            else
                logger.info(String.format("[PROGRESS] Traversed %,d %s in %.2f secs (%.2f secs per 1M %s)", nRecords, getTraversalType(), elapsed, secsPer1MReads, getTraversalType()));
        }   
    }

    /**
     * Called after a traversal to print out information about the traversal process
     */
    public void printOnTraversalDone(ReadMetrics cumulativeMetrics) {
        printProgress(null, cumulativeMetrics, true);

        final long curTime = System.currentTimeMillis();
        final double elapsed = (curTime - startTime) / 1000.0;

        // count up the number of skipped reads by summing over all filters
        long nSkippedReads = 0L;
        for ( Map.Entry<Class, Long> countsByFilter: cumulativeMetrics.getCountsByFilter().entrySet())
            nSkippedReads += countsByFilter.getValue();

        logger.info(String.format("Total runtime %.2f secs, %.2f min, %.2f hours%n", elapsed, elapsed / 60, elapsed / 3600));
        logger.info(String.format("%d reads were filtered out during traversal out of %d total (%.2f%%)",
                nSkippedReads,
                cumulativeMetrics.getNumReadsSeen(),
                100.0 * MathUtils.ratio(nSkippedReads,cumulativeMetrics.getNumReadsSeen())));
        for ( Map.Entry<Class, Long> filterCounts : cumulativeMetrics.getCountsByFilter().entrySet() ) {
            long count = filterCounts.getValue();
            logger.info(String.format("  -> %d reads (%.2f%% of total) failing %s",
                    count, 100.0 * MathUtils.ratio(count,cumulativeMetrics.getNumReadsSeen()), Utils.getClassName(filterCounts.getKey())));
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
}
