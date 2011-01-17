package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.LocusShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.ReadShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.io.DirectOutputTracker;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

import java.util.Collection;

import net.sf.picard.reference.IndexedFastaSequenceFile;


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
    public Object execute(Walker walker, ShardStrategy shardStrategy) {
        traversalEngine.startTimers();
        walker.initialize();
        Accumulator accumulator = Accumulator.create(engine,walker);

        for (Shard shard : shardStrategy) {
            if ( claimShard(shard) ) {
                // New experimental code for managing locus intervals.
                if(shard.getShardType() == Shard.ShardType.LOCUS) {
                    LocusWalker lWalker = (LocusWalker)walker;
                    WindowMaker windowMaker = new WindowMaker(shard, engine.getGenomeLocParser(), getReadIterator(shard), shard.getGenomeLocs(), lWalker.getDiscards(), engine.getSampleMetadata());
                    for(WindowMaker.WindowMakerIterator iterator: windowMaker) {
                        ShardDataProvider dataProvider = new LocusShardDataProvider(shard,iterator.getSourceInfo(),engine.getGenomeLocParser(),iterator.getLocus(),iterator,reference,rods);
                        Object result = traversalEngine.traverse(walker, dataProvider, accumulator.getReduceInit());
                        accumulator.accumulate(dataProvider,result);
                        dataProvider.close();
                    }
                    windowMaker.close();
                }
                else {
                    ShardDataProvider dataProvider = new ReadShardDataProvider(shard,engine.getGenomeLocParser(),getReadIterator(shard),reference,rods);
                    Object result = traversalEngine.traverse(walker, dataProvider, accumulator.getReduceInit());
                    accumulator.accumulate(dataProvider,result);
                    dataProvider.close();
                }
            }
        }

        Object result = accumulator.finishTraversal();

        printOnTraversalDone(result,engine.getCumulativeMetrics());

        outputTracker.close();
        cleanup();

        return accumulator;
    }

    /**
     * @{inheritDoc}
     */
    public OutputTracker getOutputTracker() { return outputTracker; }
}
