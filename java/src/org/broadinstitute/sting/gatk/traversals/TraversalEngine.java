package org.broadinstitute.sting.gatk.traversals;

import net.sf.samtools.SAMFileHeader;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;

public abstract class TraversalEngine {
    // Time in milliseconds since we initialized this engine
    private long startTime = -1;
    private long lastProgressPrintTime = -1;                // When was the last time we printed our progress?

    // How long can we go without printing some progress info?
    private final long MAX_PROGRESS_PRINT_TIME = 30 * 1000;        // 10 seconds in millisecs
    private final long N_RECORDS_TO_PRINT = 1000000;

    // Maximum number of reads to process before finishing
    protected long maxReads = -1;

    // the stored header
    private SAMFileHeader myHeader = null;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(TraversalEngine.class);

    // --------------------------------------------------------------------------------------------------------------
    //
    // Manipulating the underlying engine parameters
    //
    // --------------------------------------------------------------------------------------------------------------
    public void setMaxReads(final int maxReads) {
        this.maxReads = maxReads;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // functions for dealing locations (areas of the genome we're traversing over)
    //
    // --------------------------------------------------------------------------------------------------------------
    /**
     * get the associated SAM header for our run
     *
     * @return the header if it's stored, null if not
     */
    public SAMFileHeader getSAMHeader() {
        return myHeader;
    }

    /**
     * set's the SAM header for this traversal, which should
     * be the merged header in the multiple BAM file case.
     *
     * @param myHeader the passed in header
     */

    public void setSAMHeader(SAMFileHeader myHeader) {
        this.myHeader = myHeader;
    }
    // --------------------------------------------------------------------------------------------------------------
    //
    // printing
    //
    // --------------------------------------------------------------------------------------------------------------

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
     * @param type the type of traversal
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
     * @param <T> Type of the computation.
     */
    public abstract <T> void printOnTraversalDone(T sum);

    /**
     * Called after a traversal to print out information about the traversal process
     *
     * @param type String describing this type of traversal ("loci", "read")
     * @param sum  The reduce result of the traversal
     * @param <T>  ReduceType of the traversal
     */
    protected <T> void printOnTraversalDone(final String type, T sum) {
        printProgress(true, type, null);
        logger.info("Traversal reduce result is " + sum);
        final long curTime = System.currentTimeMillis();
        final double elapsed = (curTime - startTime) / 1000.0;
        logger.info(String.format("Total runtime %.2f secs, %.2f min, %.2f hours%n", elapsed, elapsed / 60, elapsed / 3600));
        logger.info(String.format("Traversal skipped %d valid reads out of %d total (%.2f%%)",
                TraversalStatistics.nSkippedReads,
                TraversalStatistics.nReads,
                (TraversalStatistics.nSkippedReads * 100.0) / TraversalStatistics.nReads));
        logger.info(String.format("  -> %d unmapped reads", TraversalStatistics.nUnmappedReads));
        logger.info(String.format("  -> %d duplicate reads", TraversalStatistics.nDuplicates));
        logger.info(String.format("  -> %d non-primary reads", TraversalStatistics.nNotPrimary));
        logger.info(String.format("  -> %d reads with bad alignments", TraversalStatistics.nBadAlignments));
        logger.info(String.format("  -> %d reads with indels", TraversalStatistics.nSkippedIndels));
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Initialization
    //
    // --------------------------------------------------------------------------------------------------------------

    /** Initialize the traversal engine.  After this point traversals can be run over the data */
    public void initialize() {
        lastProgressPrintTime = startTime = System.currentTimeMillis();
        // Initial the reference ordered data iterators
    }

    /**
     * this method must be implemented by all traversal engines
     *
     * @param walker       the walker to run with
     * @param shard        a shard of data
     * @param dataProvider the data provider that generates data given the shard
     * @param sum          the accumulator
     * @param <M>          an object of the map type
     * @param <T>          an object of the reduce type
     *
     * @return an object of the reduce type
     */
    public abstract <M, T> T traverse(Walker<M, T> walker,
                                      Shard shard,
                                      ShardDataProvider dataProvider,
                                      T sum);
}
