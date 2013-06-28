/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

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
        super.startingExecution();
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
