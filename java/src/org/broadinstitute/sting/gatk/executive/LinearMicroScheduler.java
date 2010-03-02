package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.io.DirectOutputTracker;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.util.Collection;

/** A micro-scheduling manager for single-threaded execution of a traversal. */
public class LinearMicroScheduler extends MicroScheduler {

    /**
     * A direct output tracker for directly managing output.
     */
    private DirectOutputTracker outputTracker = new DirectOutputTracker();

    /**
     * Create a new linear microscheduler to process the given reads and reference.
     *
     * @param walker    Walker for the traversal.
     * @param reads     Reads file(s) to process.
     * @param reference Reference for driving the traversal.
     * @param rods      Reference-ordered data.
     */
    protected LinearMicroScheduler(GenomeAnalysisEngine engine, Walker walker, SAMDataSource reads, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods ) {
        super(engine, walker, reads, reference, rods);
    }

    /**
     * Run this traversal over the specified subsection of the dataset.
     *
     * @param walker    Computation to perform over dataset.
     * @param shardStrategy A strategy for sharding the data.
     */
    public Object execute(Walker walker, ShardStrategy shardStrategy, int maxIterations) {
        // Having maxiterations in the execute method is a holdover from the old TraversalEngine days.
        // Lets do something else with this.
        traversalEngine.setMaximumIterations(maxIterations);

        walker.initialize();
        Accumulator accumulator = Accumulator.create(engine,walker);

        for (Shard shard : shardStrategy) {
            // New experimental code for managing locus intervals.
            // TODO: we'll need a similar but slightly different strategy for dealing with read intervals, so generalize this code.            
            if((shard.getShardType() == Shard.ShardType.LOCUS || shard.getShardType() == Shard.ShardType.LOCUS_INTERVAL) &&
                    shard.getGenomeLocs().size() > 0) {
                WindowMaker windowMaker = new WindowMaker(getReadIterator(shard),shard.getGenomeLocs());
                for(WindowMaker.WindowMakerIterator iterator: windowMaker) {
                    ShardDataProvider dataProvider = new ShardDataProvider(shard,iterator.getLocus(),iterator,reference,rods);
                    Object result = traversalEngine.traverse(walker, dataProvider, accumulator.getReduceInit());
                    accumulator.accumulate(dataProvider,result);
                    dataProvider.close();
                }
                windowMaker.close();
            }
            else {
                ShardDataProvider dataProvider = new ShardDataProvider(shard,null,getReadIterator(shard),reference,rods);
                Object result = traversalEngine.traverse(walker, dataProvider, accumulator.getReduceInit());
                accumulator.accumulate(dataProvider,result);
                dataProvider.close();
            }
        }

        Object result = accumulator.finishTraversal();

        printOnTraversalDone(result);

        outputTracker.close();

        return accumulator;
    }

    /**
     * @{inheritDoc}
     */
    public OutputTracker getOutputTracker() { return outputTracker; }
}
