package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;

import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;
/**
 * User: hanna
 * Date: May 18, 2009
 * Time: 2:27:17 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Manages the
 */

public abstract class Accumulator {
    /**
     * The walker for which to accumulate.
     */
    protected final Walker walker;

    /**
     * Create a new Accumulator.  Forbid outside classes from performing this operation.
     * @param walker
     */
    protected Accumulator( Walker walker ) {
        this.walker = walker;
    }

    /**
     * Creates an accumulator suitable for accumulating results of the given walker.
     * @param walker Walker for which to build an accumulator.
     * @return Accumulator suitable for this walker.s
     */
    public static Accumulator create( GenomeAnalysisEngine engine, Walker walker ) {
        if( walker.isReduceByInterval() )
            return new IntervalAccumulator( walker, engine.getIntervals() );
        else
            return new StandardAccumulator( walker );
    }

    /**
     * Gets the appropriate reduce initializer for this accumulator.
     * @return Traversal reduce init to feed into traversal engine. 
     */
    public abstract Object getReduceInit();

    /**
     * Roll this traversal result into the given accumulator.
     * @param result Result of the most recent accumulation.
     * @return the newest accumulation of the given data.
     */
    public abstract void accumulate( ShardDataProvider provider, Object result );

    /**
     * Finishes off the traversal.  Submits accumulated results to
     * the walker and returns them.
     * TODO: Its a bit funky to delegate the finishing of the traversal
     *       to an accumulator, but we're doing it for type safety so the
     *       right Walker override gets called.  Clean this up.
     * @return Final result of accumulation.
     */
    public abstract Object finishTraversal();

    /**
     * Accumulates in the 'standard' fashion; basically funnels
     * the reduce result back into the reduce init and relies on
     * the user-supplied reduce to handle the accumulation.
     */
    private static class StandardAccumulator extends Accumulator {
        private Object accumulator = null;
        private boolean initialized = false;

        protected StandardAccumulator( Walker walker ) {
            super(walker);
        }

        /**
         * Standard accumulator returns reduceInit first, then the
         * results of the previous accumulation. 
         */
        public Object getReduceInit() {
            if( !initialized ) {
                initialized = true;
                return walker.reduceInit();
            }
            else
                return accumulator;
        }

        /**
         * The result of the accumulator in a non-intervals walker
         * already takes the accumulation into account.  return the result. 
         */
        public void accumulate( ShardDataProvider provider, Object result ) { this.accumulator = result; }

        /**
         * The result of the traversal is the list of accumulated intervals.
         */
        public Object finishTraversal() {
            walker.onTraversalDone(accumulator);
            return this.accumulator;
        }
    }

    /**
     * An interval-based accumulator.  Treats each reduce result independently,
     * and aggregates those results into a single list.
     */
    private static class IntervalAccumulator extends Accumulator {
        /**
         * An iterator through all intervals in the series.
         */
        private final Iterator<GenomeLoc> intervalIterator;

        /**
         * For which interval is the accumulator currently accumulating?
         */
        private GenomeLoc currentInterval = null;

        /**
         * The actual mapping of interval to accumulator.
         */
        private final List<Pair<GenomeLoc,Object>> intervalAccumulator = new ArrayList<Pair<GenomeLoc,Object>>();

        private Object nextReduceInit = null;

        protected IntervalAccumulator(Walker walker, GenomeLocSortedSet intervals) {
            super(walker);
            this.intervalIterator = intervals.iterator();
            if(intervalIterator.hasNext()) currentInterval = intervalIterator.next();
            nextReduceInit = walker.reduceInit();
        }

        /**
         * Interval accumulator always feeds reduceInit into every new traversal.
         */
        public Object getReduceInit() { return nextReduceInit; }

        /**
         * Create a holder for interval results if none exists.  Add the result to the holder.
         */
        public void accumulate( ShardDataProvider provider, Object result ) {
            GenomeLoc location = provider.getLocus();

            // Pull the interval iterator ahead to the interval overlapping this shard fragment.
            while((currentInterval == null || currentInterval.isBefore(location)) && intervalIterator.hasNext())
                currentInterval = intervalIterator.next();

            if(currentInterval != null && currentInterval.getContig().equals(location.getContig()) && currentInterval.getStop() == location.getStop()) {
                intervalAccumulator.add(new Pair<GenomeLoc,Object>(currentInterval,result));
                nextReduceInit = walker.reduceInit();
            }
            else
                nextReduceInit = result;
        }

        /**
         * The result of the traversal is the list of accumulated intervals.
         */
        public Object finishTraversal() {
            walker.onTraversalDone(intervalAccumulator);
            return intervalAccumulator;
        }
    }
}
