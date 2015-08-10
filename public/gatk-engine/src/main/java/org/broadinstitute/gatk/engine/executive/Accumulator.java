/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.executive;

import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.datasources.providers.LocusShardDataProvider;
import org.broadinstitute.gatk.engine.datasources.providers.ShardDataProvider;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
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
        if( walker.isReduceByInterval() && engine.getIntervals() != null)
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
            walker.onTraversalDone(getReduceInit());  // must call getReduceInit to ensure that we get the accumulator value or the reduceInit value
            return this.accumulator;
        }
    }

    /**
     * An interval-based accumulator.  Treats each reduce result independently,
     * and aggregates those results into a single list.
     */
    private static class IntervalAccumulator extends Accumulator {
        /**
         * True if a new interval is being started.  This flag is used to
         * ensure that reduceInit() is not called unnecessarily.
         */
        private boolean startingNewInterval = true;

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

        /**
         * Holds the next value to be passed in as the reduce result.
         */
        private Object nextReduceInit = null;

        protected IntervalAccumulator(Walker walker, GenomeLocSortedSet intervals) {
            super(walker);
            this.intervalIterator = intervals.iterator();
            if(intervalIterator.hasNext()) currentInterval = intervalIterator.next();
        }

        /**
         * Interval accumulator always feeds reduceInit into every new traversal.
         */
        public Object getReduceInit() {
            if(startingNewInterval) {
                startingNewInterval = false;
                nextReduceInit = walker.reduceInit();
            }
            return nextReduceInit;
        }

        /**
         * Create a holder for interval results if none exists.  Add the result to the holder.
         */
        public void accumulate( ShardDataProvider provider, Object result ) {
            if(!(provider instanceof LocusShardDataProvider))
                throw new ReviewedGATKException("Unable to reduce by interval on reads traversals at this time.");

            GenomeLoc location = ((LocusShardDataProvider)provider).getLocus();

            // Pull the interval iterator ahead to the interval overlapping this shard fragment.
            while((currentInterval == null || currentInterval.isBefore(location)) && intervalIterator.hasNext())
                currentInterval = intervalIterator.next();

            if(currentInterval != null && currentInterval.getContig().equals(location.getContig()) && currentInterval.getStop() == location.getStop()) {
                intervalAccumulator.add(new Pair<GenomeLoc,Object>(currentInterval,result));
                startingNewInterval = true;
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
