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

package org.broadinstitute.gatk.engine.traversals;

import org.broadinstitute.gatk.engine.WalkerManager;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.datasources.providers.*;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.nanoScheduler.NSMapFunction;
import org.broadinstitute.gatk.utils.nanoScheduler.NSProgressFunction;
import org.broadinstitute.gatk.utils.nanoScheduler.NSReduceFunction;
import org.broadinstitute.gatk.utils.nanoScheduler.NanoScheduler;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;

import java.util.Iterator;

/**
 * A simple solution to iterating over all reference positions over a series of genomic locations.
 */
public class TraverseLociNano<M,T> extends TraversalEngine<M,T,LocusWalker<M,T>,LocusShardDataProvider> {
    /** our log, which we want to capture anything from this class */
    private static final boolean DEBUG = false;

    final NanoScheduler<MapData, MapResult, T> nanoScheduler;

    public TraverseLociNano(int nThreads) {
        nanoScheduler = new NanoScheduler<MapData, MapResult, T>(nThreads);
        nanoScheduler.setProgressFunction(new TraverseLociProgress());
    }

    @Override
    public final String getTraversalUnits() {
        return "sites";
    }

    protected static class TraverseResults<T> {
        final int numIterations;
        final T reduceResult;

        public TraverseResults(int numIterations, T reduceResult) {
            this.numIterations = numIterations;
            this.reduceResult = reduceResult;
        }
    }

    @Override
    public T traverse( LocusWalker<M,T> walker,
                       LocusShardDataProvider dataProvider,
                       T sum) {
        logger.debug(String.format("TraverseLoci.traverse: Shard is %s", dataProvider));

        final LocusView locusView = getLocusView( walker, dataProvider );

        if ( locusView.hasNext() ) { // trivial optimization to avoid unnecessary processing when there's nothing here at all
            //ReferenceOrderedView referenceOrderedDataView = new ReferenceOrderedView( dataProvider );
            ReferenceOrderedView referenceOrderedDataView = null;
            if ( WalkerManager.getWalkerDataSource(walker) != DataSource.REFERENCE_ORDERED_DATA )
                referenceOrderedDataView = new ManagingReferenceOrderedView( dataProvider );
            else
                referenceOrderedDataView = (RodLocusView)locusView;

            final LocusReferenceView referenceView = new LocusReferenceView( walker, dataProvider );

            final TraverseResults<T> result = traverse( walker, locusView, referenceView, referenceOrderedDataView, sum );
            sum = result.reduceResult;
            dataProvider.getShard().getReadMetrics().incrementNumIterations(result.numIterations);
        }

        // We have a final map call to execute here to clean up the skipped based from the
        // last position in the ROD to that in the interval
        if ( WalkerManager.getWalkerDataSource(walker) == DataSource.REFERENCE_ORDERED_DATA && ! walker.isDone() ) {
            // only do this if the walker isn't done!
            final RodLocusView rodLocusView = (RodLocusView)locusView;
            final long nSkipped = rodLocusView.getLastSkippedBases();
            if ( nSkipped > 0 ) {
                final GenomeLoc site = rodLocusView.getLocOneBeyondShard();
                final AlignmentContext ac = new AlignmentContext(site, new ReadBackedPileupImpl(site), nSkipped);
                final M x = walker.map(null, null, ac);
                sum = walker.reduce(x, sum);
            }
        }

        return sum;
    }

    /**
     * Gets the best view of loci for this walker given the available data.  The view will function as a 'trigger track'
     * of sorts, providing a consistent interface so that TraverseLoci doesn't need to be reimplemented for any new datatype
     * that comes along.
     * @param walker walker to interrogate.
     * @param dataProvider Data which which to drive the locus view.
     * @return A view of the locus data, where one iteration of the locus view maps to one iteration of the traversal.
     */
    private LocusView getLocusView( Walker<M,T> walker, LocusShardDataProvider dataProvider ) {
        final DataSource dataSource = WalkerManager.getWalkerDataSource(walker);
        if( dataSource == DataSource.READS )
            return new CoveredLocusView(dataProvider);
        else if( dataSource == DataSource.REFERENCE ) //|| ! GenomeAnalysisEngine.instance.getArguments().enableRodWalkers )
            return new AllLocusView(dataProvider);
        else if( dataSource == DataSource.REFERENCE_ORDERED_DATA )
            return new RodLocusView(dataProvider);
        else
            throw new UnsupportedOperationException("Unsupported traversal type: " + dataSource);
    }

    protected TraverseResults<T> traverse(final LocusWalker<M, T> walker,
                                          final LocusView locusView,
                                          final LocusReferenceView referenceView,
                                          final ReferenceOrderedView referenceOrderedDataView,
                                          final T sum) {
        nanoScheduler.setDebug(DEBUG);
        final TraverseLociMap myMap = new TraverseLociMap(walker);
        final TraverseLociReduce myReduce = new TraverseLociReduce(walker);

        final MapDataIterator inputIterator = new MapDataIterator(locusView, referenceView, referenceOrderedDataView);
        final T result = nanoScheduler.execute(inputIterator, myMap, sum, myReduce);

        return new TraverseResults<T>(inputIterator.numIterations, result);
    }

    /**
     * Create iterator that provides inputs for all map calls into MapData, to be provided
     * to NanoScheduler for Map/Reduce
     */
    private class MapDataIterator implements Iterator<MapData> {
        final LocusView locusView;
        final LocusReferenceView referenceView;
        final ReferenceOrderedView referenceOrderedDataView;
        int numIterations = 0;

        private MapDataIterator(LocusView locusView, LocusReferenceView referenceView, ReferenceOrderedView referenceOrderedDataView) {
            this.locusView = locusView;
            this.referenceView = referenceView;
            this.referenceOrderedDataView = referenceOrderedDataView;
        }

        @Override
        public boolean hasNext() {
            return locusView.hasNext() && ! engine.exceedsRuntimeLimit();
        }

        @Override
        public MapData next() {
            final AlignmentContext locus = locusView.next();
            final GenomeLoc location = locus.getLocation();

            //logger.info("Pulling data from MapDataIterator at " + location);

            // create reference context. Note that if we have a pileup of "extended events", the context will
            // hold the (longest) stretch of deleted reference bases (if deletions are present in the pileup).
            final ReferenceContext refContext = referenceView.getReferenceContext(location);

            // Iterate forward to get all reference ordered data covering this location
            final RefMetaDataTracker tracker = referenceOrderedDataView.getReferenceOrderedDataAtLocus(location);

            numIterations++;
            return new MapData(locus, refContext,  tracker);
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Cannot remove elements from MapDataIterator");
        }
    }

    @Override
    public void shutdown() {
        nanoScheduler.shutdown();
    }

    /**
     * The input data needed for each map call.  The read, the reference, and the RODs
     */
    private class MapData {
        final AlignmentContext alignmentContext;
        final ReferenceContext refContext;
        final RefMetaDataTracker tracker;

        private MapData(final AlignmentContext alignmentContext, ReferenceContext refContext, RefMetaDataTracker tracker) {
            this.alignmentContext = alignmentContext;
            this.refContext = refContext;
            this.tracker = tracker;
        }

        @Override
        public String toString() {
            return "MapData " + alignmentContext.getLocation();
        }
    }

    /**
     * Contains the results of a map call, indicating whether the call was good, filtered, or done
     */
    private class MapResult {
        final M value;
        final boolean reduceMe;

        /**
         * Create a MapResult with value that should be reduced
         *
         * @param value the value to reduce
         */
        private MapResult(final M value) {
            this.value = value;
            this.reduceMe = true;
        }

        /**
         * Create a MapResult that shouldn't be reduced
         */
        private MapResult() {
            this.value = null;
            this.reduceMe = false;
        }
    }

    /**
     * A static object that tells reduce that the result of map should be skipped (filtered or done)
     */
    private final MapResult SKIP_REDUCE = new MapResult();

    /**
     * MapFunction for TraverseReads meeting NanoScheduler interface requirements
     *
     * Applies walker.map to MapData, returning a MapResult object containing the result
     */
    private class TraverseLociMap implements NSMapFunction<MapData, MapResult> {
        final LocusWalker<M,T> walker;

        private TraverseLociMap(LocusWalker<M, T> walker) {
            this.walker = walker;
        }

        @Override
        public MapResult apply(final MapData data) {
            if ( ! walker.isDone() ) {
                final boolean keepMeP = walker.filter(data.tracker, data.refContext, data.alignmentContext);
                if (keepMeP) {
                    final M x = walker.map(data.tracker, data.refContext, data.alignmentContext);
                    return new MapResult(x);
                }
            }
            return SKIP_REDUCE;
        }
    }

    /**
     * NSReduceFunction for TraverseReads meeting NanoScheduler interface requirements
     *
     * Takes a MapResult object and applies the walkers reduce function to each map result, when applicable
     */
    private class TraverseLociReduce implements NSReduceFunction<MapResult, T> {
        final LocusWalker<M,T> walker;

        private TraverseLociReduce(LocusWalker<M, T> walker) {
            this.walker = walker;
        }

        @Override
        public T apply(MapResult one, T sum) {
            if ( one.reduceMe )
                // only run reduce on values that aren't DONE or FAILED
                return walker.reduce(one.value, sum);
            else
                return sum;
        }
    }

    private class TraverseLociProgress implements NSProgressFunction<MapData> {
        @Override
        public void progress(MapData lastProcessedMap) {
            if (lastProcessedMap.alignmentContext != null)
                printProgress(lastProcessedMap.alignmentContext.getLocation());
        }
    }
}
