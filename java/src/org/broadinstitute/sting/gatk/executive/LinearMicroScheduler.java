package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.dataSources.providers.LocusContextProvider;
import org.broadinstitute.sting.gatk.dataSources.providers.ReferenceProvider;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SimpleDataSourceLoadException;
import org.broadinstitute.sting.gatk.iterators.MergingSamRecordIterator2;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.gatk.traversals.TraverseLociByReference;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.File;
import java.util.List;

/**
 * A micro-scheduling manager for single-threaded execution of a traversal.
 */
public class LinearMicroScheduler extends MicroScheduler {

    private TraverseLociByReference traversalEngine = null;

    public TraversalEngine getTraversalEngine() {
        return traversalEngine;
    }

    /**
     * Create a new linear microscheduler to process the given reads and reference.
     * @param reads Reads file(s) to process.
     * @param refFile Reference for driving the traversal.
     */
    protected LinearMicroScheduler( List<File> reads, File refFile ) {
        super( reads, refFile );
        traversalEngine = new TraverseLociByReference( reads, refFile, new java.util.ArrayList() );
    }

    /**
     * Run this traversal over the specified subsection of the dataset.
     * @param walker Computation to perform over dataset.
     * @param locations Subset of the dataset over which to walk.
     */
    public void execute( Walker walker, List<GenomeLoc> locations ) {
        ShardStrategy shardStrategy = getShardStrategy( reference, locations );
        SAMDataSource dataSource = getReadsDataSource();

        boolean walkerInitialized = false;
        Object accumulator = null;

        for(Shard shard: shardStrategy) {
            GenomeLoc span = shard.getGenomeLoc();

            MergingSamRecordIterator2 readShard = null;
            try {
                readShard = dataSource.seek( span );
            }
            catch( SimpleDataSourceLoadException ex ) {
                throw new RuntimeException( ex );
            }

            ReferenceProvider referenceProvider = new ReferenceProvider( reference, span );
            LocusContextProvider locusProvider = new LocusContextProvider( readShard );

            // set the sam header of the traversal engine
            traversalEngine.setSAMHeader(readShard.getMergedHeader());

            if (!walkerInitialized) {
                walker.initialize();
                accumulator = ((LocusWalker<?,?>)walker).reduceInit();
                walkerInitialized = true;
            }

            accumulator = traversalEngine.traverse( walker, shard, referenceProvider, locusProvider, accumulator );
            readShard.close();
        }

        traversalEngine.printOnTraversalDone("loci", accumulator);
        walker.onTraversalDone(accumulator);
    }


}
