package org.broadinstitute.sting.gatk.executive;

import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.dataSources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.gatk.traversals.TraverseReads;
import org.broadinstitute.sting.gatk.traversals.TraverseLociByReference;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
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
    private static long SHARD_SIZE = 100000L;    

    protected static Logger logger = Logger.getLogger(MicroScheduler.class);

    protected final TraversalEngine traversalEngine;
    protected final IndexedFastaSequenceFile reference;

    private final SAMDataSource reads;

    /**
     * MicroScheduler factory function.  Create a microscheduler appropriate for reducing the
     * selected walker.
     * @param walker Which walker to use.
     * @param nThreadsToUse Number of threads to utilize.
     * @return The best-fit microscheduler.
     */
    public static MicroScheduler create( Walker walker, Reads reads, File ref, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods, int nThreadsToUse ) {
        if( walker instanceof TreeReducible && nThreadsToUse > 1 ) {
            logger.info("Creating hierarchical microscheduler");
            return new HierarchicalMicroScheduler( walker, reads, ref, rods, nThreadsToUse );
        }
        else {
            logger.info("Creating linear microscheduler");
            return new LinearMicroScheduler( walker, reads, ref, rods );
        }
    }

    /**
     * Create a microscheduler given the reads and reference.
     * @param reads The reads.
     * @param refFile File pointer to the reference.
     */
    protected MicroScheduler( Walker walker, Reads reads, File refFile, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods ) {
        if (walker instanceof ReadWalker) {
            traversalEngine = new TraverseReads(reads.getReadsFiles(), refFile, rods);
        } else {
            traversalEngine = new TraverseLociByReference(reads.getReadsFiles(), refFile, rods);
        }

        this.reads = getReadsDataSource( reads );
        this.reference = openReferenceSequenceFile( refFile );
    }

    /**
     * A temporary getter for the traversal engine.  In the future, clients
     * of the microscheduler shouldn't need to know anything about the traversal engine.
     * @return The traversal engine.
     */
    public TraversalEngine getTraversalEngine() {
        return traversalEngine;
    }

    /**
     * Walks a walker over the given list of intervals.
     * @param walker Computation to perform over dataset.
     * @param intervals A list of intervals over which to walk.  Null for whole dataset.
     * @return the return type of the walker
     */
    public abstract Object execute( Walker walker, List<GenomeLoc> intervals);

    /**
     * Get the sharding strategy given a driving data source.
     * @param drivingDataSource Data on which to shard.
     * @param intervals Intervals to use when limiting sharding.
     * @return Sharding strategy for this driving data source.
     */
    protected ShardStrategy getShardStrategy( ReferenceSequenceFile drivingDataSource, List<GenomeLoc> intervals ) {
        ShardStrategy shardStrategy = null;
        ShardStrategyFactory.SHATTER_STRATEGY shardType = (traversalEngine instanceof TraverseReads) ?
                ShardStrategyFactory.SHATTER_STRATEGY.READS : ShardStrategyFactory.SHATTER_STRATEGY.LINEAR;
        
        if( intervals != null && shardType != ShardStrategyFactory.SHATTER_STRATEGY.READS)
            shardStrategy = ShardStrategyFactory.shatter( shardType,
                                                          drivingDataSource.getSequenceDictionary(),
                                                          SHARD_SIZE,
                                                          intervals );
        else
            shardStrategy = ShardStrategyFactory.shatter( shardType,
                                                          drivingDataSource.getSequenceDictionary(),
                                                          SHARD_SIZE );

        return shardStrategy;
    }

    /**
     * Gets an window into all the data that can be viewed as a single shard.
     * @param shard The section of data to view.
     * @return An accessor for all the data in this shard.
     */
    protected ShardDataProvider getShardDataProvider( Shard shard ) {
        return new ShardDataProvider( shard, reads, reference );
    }

    /**
     * Gets a data source for the given set of reads.
     * @return A data source for the given set of reads.
     */
    private SAMDataSource getReadsDataSource( Reads reads ) {
        SAMDataSource dataSource = new SAMDataSource( reads );

        // Side effect: initialize the traversal engine with reads data.
        // TODO: Give users a dedicated way of getting the header so that the MicroScheduler
        //       doesn't have to bend over backward providing legacy getters and setters.
        traversalEngine.setSAMHeader(dataSource.getHeader());

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
