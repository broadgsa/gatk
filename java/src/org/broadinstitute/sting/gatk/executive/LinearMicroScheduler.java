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
 * A micro-scheduling manager for N-way threaded execution of a traversal
 *
 */
public class LinearMicroScheduler extends MicroScheduler {

    private TraverseLociByReference traversalEngine = null;

    public TraversalEngine getTraversalEngine() {
        return traversalEngine;
    }

    public LinearMicroScheduler( List<File> reads,              // the reads file(s)
                         File refFile,                        // the reference file driving the traversal
                         int nThreadsToUse ) {                // maximum number of threads to use to do the work
        super( reads, refFile );
        traversalEngine = new TraverseLociByReference( reads, refFile, new java.util.ArrayList() );
    }

    public void execute( Walker walker,                        // the analysis technique to use.
                         List<GenomeLoc> locations ) {         // list of work to do

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
