package org.broadinstitute.sting.gatk.executive;

import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SimpleDataSourceLoadException;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 26, 2009
 * Time: 12:37:23 PM
 * To change this template use File | Settings | File Templates.
 */

/**
 * Shards and schedules data in manageable chunks.
 */
public abstract class MicroScheduler {
    private List<File> reads;
    private static long SHARD_SIZE = 100000L;    

    protected static Logger logger = Logger.getLogger(MicroScheduler.class);

    protected IndexedFastaSequenceFile reference;

    /**
     * MicroScheduler factory function.  Create a microscheduler appropriate for reducing the
     * selected walker.
     * @param walker Which walker to use.
     * @param nThreadsToUse Number of threads to utilize.
     * @return The best-fit microscheduler.
     */
    public static MicroScheduler create( Walker walker, List<File> reads, File ref, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods, int nThreadsToUse ) {
        if( walker instanceof TreeReducible && nThreadsToUse > 1 ) {
            logger.info("Creating hierarchical microscheduler");
            return new HierarchicalMicroScheduler( reads, ref, rods, nThreadsToUse );
        }
        else {
            logger.info("Creating linear microscheduler");
            return new LinearMicroScheduler( reads, ref, rods, walker );
        }
    }

    /**
     * Create a microscheduler given the reads and reference.
     * @param reads The reads.
     * @param refFile File pointer to the reference.
     */
    protected MicroScheduler( List<File> reads, File refFile ) {
        this.reads = reads;
        this.reference = openReferenceSequenceFile( refFile );
    }

    /**
     * A temporary getter for the traversal engine.  In the future, clients
     * of the microscheduler shouldn't need to know anything about the traversal engine.
     * @return The traversal engine.
     */
    public abstract TraversalEngine getTraversalEngine();

    /**
     * Walks a walker over the given list of intervals.
     * @param walker Computation to perform over dataset.
     * @param intervals A list of intervals over which to walk.  Null for whole dataset.
     */
    public abstract void execute( Walker walker, List<GenomeLoc> intervals);

    /**
     * Get the sharding strategy given a driving data source.
     * @param drivingDataSource Data on which to shard.
     * @param intervals Intervals to use when limiting sharding.
     * @return Sharding strategy for this driving data source.
     */
    protected ShardStrategy getShardStrategy( ReferenceSequenceFile drivingDataSource, List<GenomeLoc> intervals ) {
        ShardStrategy shardStrategy = null;
        if( intervals != null )
            shardStrategy = ShardStrategyFactory.shatter( ShardStrategyFactory.SHATTER_STRATEGY.LINEAR,
                                                          drivingDataSource.getSequenceDictionary(),
                                                          SHARD_SIZE,
                                                          intervals );
        else
            shardStrategy = ShardStrategyFactory.shatter( ShardStrategyFactory.SHATTER_STRATEGY.LINEAR,
                                                          drivingDataSource.getSequenceDictionary(),
                                                          SHARD_SIZE );

        return shardStrategy;
    }

    /**
     * Gets a data source for the given set of reads.
     * @return A data source for the given set of reads.
     */
    protected SAMDataSource getReadsDataSource() {
        SAMDataSource dataSource = null;
        try {
            dataSource = new SAMDataSource( TraversalEngine.unpackReads(reads) );
        }
        catch( SimpleDataSourceLoadException ex ) {
            throw new RuntimeException( ex );
        }
        catch( FileNotFoundException ex ) {
            throw new RuntimeException( ex );
        }
        return dataSource;
    }

    /**
     * Opens a reference sequence file paired with an index.
     * @param refFile Handle to a reference sequence file.  Non-null.
     * @return A thread-safe file wrapper.
     */
    private IndexedFastaSequenceFile openReferenceSequenceFile( File refFile ) {
        IndexedFastaSequenceFile ref = null;
        try {
            ref = new IndexedFastaSequenceFile(refFile);
        }
        catch( FileNotFoundException ex ) {
            throw new RuntimeException("File not found opening fasta file; please do this check before MicroManaging", ex);
        }
        GenomeLoc.setupRefContigOrdering(ref);
        return ref;
    }
}
