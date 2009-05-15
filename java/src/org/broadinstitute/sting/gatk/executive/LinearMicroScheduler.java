package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.dataSources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.File;
import java.util.List;

/** A micro-scheduling manager for single-threaded execution of a traversal. */
public class LinearMicroScheduler extends MicroScheduler {

    /**
     * Create a new linear microscheduler to process the given reads and reference.
     *
     * @param reads   Reads file(s) to process.
     * @param refFile Reference for driving the traversal.
     */
    protected LinearMicroScheduler( Walker walker, List<File> reads, File refFile, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods ) {
        super(walker, reads, refFile, rods);
    }

    /**
     * Run this traversal over the specified subsection of the dataset.
     *
     * @param walker    Computation to perform over dataset.
     * @param locations Subset of the dataset over which to walk.
     */
    public Object execute(Walker walker, List<GenomeLoc> locations) {
        ShardStrategy shardStrategy = getShardStrategy(reference, locations);

        walker.initialize();
        Object accumulator = walker.reduceInit();

        for (Shard shard : shardStrategy) {
            ShardDataProvider dataProvider = getShardDataProvider( shard );
            accumulator = traversalEngine.traverse(walker, shard, dataProvider, accumulator);
            dataProvider.close();
        }

        traversalEngine.printOnTraversalDone(accumulator);

        walker.onTraversalDone(accumulator);

        return accumulator;
    }


}
