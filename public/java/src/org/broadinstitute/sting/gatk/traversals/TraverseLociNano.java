package org.broadinstitute.sting.gatk.traversals;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.LocusReferenceView;
import org.broadinstitute.sting.gatk.datasources.providers.LocusView;
import org.broadinstitute.sting.gatk.datasources.providers.ReferenceOrderedView;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.nanoScheduler.MapFunction;
import org.broadinstitute.sting.utils.nanoScheduler.NanoScheduler;
import org.broadinstitute.sting.utils.nanoScheduler.ReduceFunction;

import java.util.Iterator;

/**
 * A simple solution to iterating over all reference positions over a series of genomic locations.
 */
public class TraverseLociNano<M,T> extends TraverseLociBase<M,T> {
    /** our log, which we want to capture anything from this class */
    private static final boolean DEBUG = false;
    private static final int BUFFER_SIZE = 1000;

    final NanoScheduler<MapData, MapResult, T> nanoScheduler;

    public TraverseLociNano(int nThreads) {
        nanoScheduler = new NanoScheduler<MapData, MapResult, T>(BUFFER_SIZE, nThreads);
    }

    @Override
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

        // todo -- how do I print progress?
//        final GATKSAMRecord lastRead = aggregatedInputs.get(aggregatedInputs.size() - 1).read;
//        final GenomeLoc locus = engine.getGenomeLocParser().createGenomeLoc(lastRead);
//        printProgress(dataProvider.getShard(), locus, aggregatedInputs.size());

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
            return locusView.hasNext();
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
            final RefMetaDataTracker tracker = referenceOrderedDataView.getReferenceOrderedDataAtLocus(location, refContext);

            numIterations++;
            return new MapData(locus, refContext,  tracker);
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Cannot remove elements from MapDataIterator");
        }
    }

    @Override
    public void printOnTraversalDone() {
        nanoScheduler.shutdown();
        super.printOnTraversalDone();
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
    private class TraverseLociMap implements MapFunction<MapData, MapResult> {
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
     * ReduceFunction for TraverseReads meeting NanoScheduler interface requirements
     *
     * Takes a MapResult object and applies the walkers reduce function to each map result, when applicable
     */
    private class TraverseLociReduce implements ReduceFunction<MapResult, T> {
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
}
