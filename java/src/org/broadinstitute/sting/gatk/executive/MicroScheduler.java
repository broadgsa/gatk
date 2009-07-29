/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.executive;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.traversals.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.util.*;


/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 26, 2009
 * Time: 12:37:23 PM
 * To change this template use File | Settings | File Templates.
 */

/** Shards and schedules data in manageable chunks. */
public abstract class MicroScheduler {
    protected static Logger logger = Logger.getLogger(MicroScheduler.class);

    protected final TraversalEngine traversalEngine;
    protected final IndexedFastaSequenceFile reference;

    private final SAMDataSource reads;
    private final Collection<ReferenceOrderedDataSource> rods;

    /**
     * MicroScheduler factory function.  Create a microscheduler appropriate for reducing the
     * selected walker.
     *
     * @param walker        Which walker to use.
     * @param reads         the informations associated with the reads
     * @param reference     the reference file
     * @param rods          the rods to include in the traversal
     * @param nThreadsToUse Number of threads to utilize.
     *
     * @return The best-fit microscheduler.
     */
    public static MicroScheduler create(Walker walker, SAMDataSource reads, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods, int nThreadsToUse) {
        if (walker instanceof TreeReducible && nThreadsToUse > 1) {
            logger.info("Creating hierarchical microscheduler");
            return new HierarchicalMicroScheduler(walker, reads, reference, rods, nThreadsToUse);
        } else {
            logger.info("Creating linear microscheduler");
            return new LinearMicroScheduler(walker, reads, reference, rods);
        }
    }

    /**
     * Create a microscheduler given the reads and reference.
     *
     * @param walker  the walker to execute with
     * @param reads   The reads.
     * @param reference The reference.
     * @param rods    the rods to include in the traversal
     */
    protected MicroScheduler(Walker walker, SAMDataSource reads, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods) {
        this.reads = reads;
        this.reference = reference;
        this.rods = rods;

        if (walker instanceof ReadWalker) {
            traversalEngine = new TraverseReads();
        } else if (walker instanceof LocusWalker) {
            traversalEngine = new TraverseLoci();
        } else if (walker instanceof LocusWindowWalker) {
            traversalEngine = new TraverseLocusWindows();
        } else if (walker instanceof DuplicateWalker) {
            traversalEngine = new TraverseDuplicates();
        } else {
            throw new UnsupportedOperationException("Unable to determine traversal type, the walker is an unknown type.");
        }        

        traversalEngine.initialize();                
    }

    /**
     * Walks a walker over the given list of intervals.
     *
     * @param walker        Computation to perform over dataset.
     * @param shardStrategy A strategy for sharding the data.
     *
     * @return the return type of the walker
     */
    public abstract Object execute(Walker walker, ShardStrategy shardStrategy, int iterations );


    /**
     * Gets an window into all the data that can be viewed as a single shard.
     *
     * @param shard The section of data to view.
     *
     * @return An accessor for all the data in this shard.
     */
    protected ShardDataProvider getShardDataProvider(Shard shard) {
        return new ShardDataProvider(shard, reads, reference, rods);
    }

    /**
     * Print summary information for the analysis.
     * @param sum The final reduce output.
     */
    protected void printOnTraversalDone(Object sum) {
        // HACK: The microscheduler should be too dumb to know anything about the data
        //       it's actually processing; it should just funnel anything it receives
        //       to the traversal engine.
        // TODO: Implement code to allow the datasources to print summary info of the
        //       data they've seen.
        if( reads != null && reads.getViolationHistogram().getViolationCount() > 0 )
            logger.warn(String.format("%n%s",reads.getViolationHistogram()));

        traversalEngine.printOnTraversalDone(sum);
    }

    /**
     * Returns data source maintained by this scheduler
     * @return
     */
    public SAMDataSource getSAMDataSource() { return reads; }

    /**
     * Returns the reference maintained by this scheduler.
     * @return The reference maintained by this scheduler.
     */
    public IndexedFastaSequenceFile getReference() { return reference; }
}
