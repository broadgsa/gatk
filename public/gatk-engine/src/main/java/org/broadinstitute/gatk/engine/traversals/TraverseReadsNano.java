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

import htsjdk.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.datasources.providers.ReadBasedReferenceOrderedView;
import org.broadinstitute.gatk.engine.datasources.providers.ReadReferenceView;
import org.broadinstitute.gatk.engine.datasources.providers.ReadShardDataProvider;
import org.broadinstitute.gatk.engine.datasources.providers.ReadView;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.utils.nanoScheduler.NSMapFunction;
import org.broadinstitute.gatk.utils.nanoScheduler.NSProgressFunction;
import org.broadinstitute.gatk.utils.nanoScheduler.NSReduceFunction;
import org.broadinstitute.gatk.utils.nanoScheduler.NanoScheduler;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.Iterator;
import java.util.LinkedList;

/**
 * A nano-scheduling version of TraverseReads.
 *
 * Implements the traversal of a walker that accepts individual reads, the reference, and
 * RODs per map call.  Directly supports shared memory parallelism via NanoScheduler
 *
 * @author depristo
 * @version 1.0
 * @date 9/2/2012
 */
public class TraverseReadsNano<M,T> extends TraversalEngine<M,T,ReadWalker<M,T>,ReadShardDataProvider> {
    /** our log, which we want to capture anything from this class */
    private final static boolean PRE_READ_ALL_MAP_DATA = true;
    protected static final Logger logger = Logger.getLogger(TraverseReadsNano.class);
    private static final boolean DEBUG = false;
    final NanoScheduler<MapData, MapResult, T> nanoScheduler;

    public TraverseReadsNano(int nThreads) {
        nanoScheduler = new NanoScheduler<MapData, MapResult, T>(nThreads);
        nanoScheduler.setProgressFunction(new NSProgressFunction<MapData>() {
            @Override
            public void progress(MapData lastProcessedMap) {
                if ( lastProcessedMap.refContext != null )
                    // note, need to use getStopLocation so we don't give an interval to ProgressMeterDaemon
                    printProgress(lastProcessedMap.refContext.getLocus().getStopLocation());
            }
        });
    }

    @Override
    public String getTraversalUnits() {
        return "reads";
    }

    /**
     * Traverse by reads, given the data and the walker
     *
     * @param walker the walker to traverse with
     * @param dataProvider the provider of the reads data
     * @param sum the value of type T, specified by the walker, to feed to the walkers reduce function
     * @return the reduce variable of the read walker
     */
    public T traverse(ReadWalker<M,T> walker,
                      ReadShardDataProvider dataProvider,
                      T sum) {
        if ( logger.isDebugEnabled() )
            logger.debug(String.format("TraverseReadsNano.traverse Covered dataset is %s", dataProvider));

        if( !dataProvider.hasReads() )
            throw new IllegalArgumentException("Unable to traverse reads; no read data is available.");

        nanoScheduler.setDebug(DEBUG);
        final TraverseReadsMap myMap = new TraverseReadsMap(walker);
        final TraverseReadsReduce myReduce = new TraverseReadsReduce(walker);

        final Iterator<MapData> aggregatedInputs = aggregateMapData(dataProvider);
        final T result = nanoScheduler.execute(aggregatedInputs, myMap, sum, myReduce);

        return result;
    }

    /**
     * Aggregate all of the inputs for all map calls into MapData, to be provided
     * to NanoScheduler for Map/Reduce
     *
     * @param dataProvider the source of our data
     * @return a linked list of MapData objects holding the read, ref, and ROD info for every map/reduce
     *          should execute
     */
    private Iterator<MapData> aggregateMapData(final ReadShardDataProvider dataProvider) {
        final Iterator<MapData> it = makeDataIterator(dataProvider);
        if ( PRE_READ_ALL_MAP_DATA ) {
            final LinkedList<MapData> l = new LinkedList<MapData>();
            while ( it.hasNext() ) l.add(it.next());
            return l.iterator();
        } else {
            return it;
        }
    }


    private Iterator<MapData> makeDataIterator(final ReadShardDataProvider dataProvider) {
        return new Iterator<MapData> ()  {
            final ReadView reads = new ReadView(dataProvider);
            final ReadReferenceView reference = new ReadReferenceView(dataProvider);
            final ReadBasedReferenceOrderedView rodView = new ReadBasedReferenceOrderedView(dataProvider);
            final Iterator<SAMRecord> readIterator = reads.iterator();

            @Override public boolean hasNext() { return ! engine.exceedsRuntimeLimit() && readIterator.hasNext(); }

            @Override
            public MapData next() {
                final SAMRecord read = readIterator.next();
                final ReferenceContext refContext = ! read.getReadUnmappedFlag()
                        ? reference.getReferenceContext(read)
                        : null;

                // if the read is mapped, create a metadata tracker
                final RefMetaDataTracker tracker = read.getReferenceIndex() >= 0
                        ? rodView.getReferenceOrderedDataForRead(read)
                        : null;

                // update the number of reads we've seen
                dataProvider.getShard().getReadMetrics().incrementNumIterations();

                return new MapData((GATKSAMRecord)read, refContext, tracker);
            }

            @Override public void remove() {
                throw new UnsupportedOperationException("Remove not supported");
            }
        };
    }

    @Override
    public void shutdown() {
        nanoScheduler.shutdown();
    }

    /**
     * The input data needed for each map call.  The read, the reference, and the RODs
     */
    private class MapData {
        final GATKSAMRecord read;
        final ReferenceContext refContext;
        final RefMetaDataTracker tracker;

        private MapData(GATKSAMRecord read, ReferenceContext refContext, RefMetaDataTracker tracker) {
            this.read = read;
            this.refContext = refContext;
            this.tracker = tracker;
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
    private class TraverseReadsMap implements NSMapFunction<MapData, MapResult> {
        final ReadWalker<M,T> walker;

        private TraverseReadsMap(ReadWalker<M, T> walker) {
            this.walker = walker;
        }

        @Override
        public MapResult apply(final MapData data) {
            if ( ! walker.isDone() ) {
                final boolean keepMeP = walker.filter(data.refContext, data.read);
                if (keepMeP)
                    return new MapResult(walker.map(data.refContext, data.read, data.tracker));
            }

            return SKIP_REDUCE;
        }
    }

    /**
     * NSReduceFunction for TraverseReads meeting NanoScheduler interface requirements
     *
     * Takes a MapResult object and applies the walkers reduce function to each map result, when applicable
     */
    private class TraverseReadsReduce implements NSReduceFunction<MapResult, T> {
        final ReadWalker<M,T> walker;

        private TraverseReadsReduce(ReadWalker<M, T> walker) {
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
