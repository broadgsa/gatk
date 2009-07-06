package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;

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
    protected LinearMicroScheduler( Walker walker, Reads reads, File refFile, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods ) {
        super(walker, reads, refFile, rods);
    }

    /**
     * Run this traversal over the specified subsection of the dataset.
     *
     * @param walker    Computation to perform over dataset.
     * @param locations Subset of the dataset over which to walk.
     * @param maxIterations the maximum number of iterations we're to perform
     */
    public Object execute(Walker walker, GenomeLocSortedSet locations, Integer maxIterations) {
        ShardStrategy shardStrategy = getShardStrategy(walker, reference, locations, maxIterations);

        walker.initialize();
        Accumulator accumulator = Accumulator.create(walker);

        for (Shard shard : shardStrategy) {
            ShardDataProvider dataProvider = getShardDataProvider( shard );

            Object result = traversalEngine.traverse(walker, shard, dataProvider, accumulator.getReduceInit());
            accumulator.accumulate( shard, result );

            dataProvider.close();
        }

        Object result = accumulator.finishTraversal();

        printOnTraversalDone(result);

        return accumulator;
    }


}
