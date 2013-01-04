package org.broadinstitute.sting.gatk.executive;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.datasources.providers.LocusShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.ReadShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.reads.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.io.DirectOutputTracker;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.gatk.resourcemanagement.ThreadAllocation;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.gatk.traversals.TraverseActiveRegions;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.threading.ThreadEfficiencyMonitor;

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
    protected LinearMicroScheduler(final GenomeAnalysisEngine engine,
                                   final Walker walker,
                                   final SAMDataSource reads,
                                   final IndexedFastaSequenceFile reference,
                                   final Collection<ReferenceOrderedDataSource> rods,
                                   final ThreadAllocation threadAllocation) {
        super(engine, walker, reads, reference, rods, threadAllocation);

        if ( threadAllocation.monitorThreadEfficiency() )
            setThreadEfficiencyMonitor(new ThreadEfficiencyMonitor());
    }

    /**
     * Run this traversal over the specified subsection of the dataset.
     *
     * @param walker    Computation to perform over dataset.
     * @param shardStrategy A strategy for sharding the data.
     */
    public Object execute(Walker walker, Iterable<Shard> shardStrategy) {
        walker.initialize();
        Accumulator accumulator = Accumulator.create(engine,walker);

        boolean done = walker.isDone();
        int counter = 0;

        final TraversalEngine traversalEngine = borrowTraversalEngine(this);
        for (Shard shard : shardStrategy ) {
            if ( abortExecution() || done || shard == null ) // we ran out of shards that aren't owned
                break;

            if(shard.getShardType() == Shard.ShardType.LOCUS) {
                WindowMaker windowMaker = new WindowMaker(shard, engine.getGenomeLocParser(),
                        getReadIterator(shard), shard.getGenomeLocs(), SampleUtils.getSAMFileSamples(engine));
                for(WindowMaker.WindowMakerIterator iterator: windowMaker) {
                    ShardDataProvider dataProvider = new LocusShardDataProvider(shard,iterator.getSourceInfo(),engine.getGenomeLocParser(),iterator.getLocus(),iterator,reference,rods);
                    Object result = traversalEngine.traverse(walker, dataProvider, accumulator.getReduceInit());
                    accumulator.accumulate(dataProvider,result);
                    dataProvider.close();
                    if ( walker.isDone() ) break;
                }
                windowMaker.close();
            }
            else {
                ShardDataProvider dataProvider = new ReadShardDataProvider(shard,engine.getGenomeLocParser(),getReadIterator(shard),reference,rods);
                Object result = traversalEngine.traverse(walker, dataProvider, accumulator.getReduceInit());
                accumulator.accumulate(dataProvider,result);
                dataProvider.close();
            }

            done = walker.isDone();
        }

        // Special function call to empty out the work queue. Ugly for now but will be cleaned up when we eventually push this functionality more into the engine
        if( traversalEngine instanceof TraverseActiveRegions ) {
            final Object result = ((TraverseActiveRegions) traversalEngine).endTraversal(walker, accumulator.getReduceInit());
            accumulator.accumulate(null, result); // Assumes only used with StandardAccumulator
        }
                
        Object result = accumulator.finishTraversal();

        outputTracker.close();
        returnTraversalEngine(this, traversalEngine);
        cleanup();
        executionIsDone();

        return accumulator;
    }

    /**
     * @{inheritDoc}
     */
    public OutputTracker getOutputTracker() { return outputTracker; }
}
