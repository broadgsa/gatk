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
package org.broadinstitute.sting.gatk.traversals;

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.ReadBasedReferenceOrderedView;
import org.broadinstitute.sting.gatk.datasources.providers.ReadReferenceView;
import org.broadinstitute.sting.gatk.datasources.providers.ReadShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.ReadView;
import org.broadinstitute.sting.gatk.datasources.reads.ReadShard;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.nanoScheduler.MapFunction;
import org.broadinstitute.sting.utils.nanoScheduler.NanoScheduler;
import org.broadinstitute.sting.utils.nanoScheduler.ReduceFunction;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.ArrayList;
import java.util.List;

/**
 * @author aaron
 * @version 1.0
 * @date Apr 24, 2009
 * <p/>
 * Class TraverseReads
 * <p/>
 * This class handles traversing by reads in the new shardable style
 */
public class TraverseReadsNano<M,T> extends TraversalEngine<M,T,ReadWalker<M,T>,ReadShardDataProvider> {
    /** our log, which we want to capture anything from this class */
    protected static final Logger logger = Logger.getLogger(TraverseReadsNano.class);
    private static final boolean DEBUG = false;
    final NanoScheduler<MapData, M, T> nanoScheduler;

    public TraverseReadsNano(int nThreads) {
        final int bufferSize = ReadShard.getReadBufferSize() + 1; // actually has 1 more than max
        nanoScheduler = new NanoScheduler<MapData, M, T>(bufferSize, nThreads);
    }

    @Override
    protected String getTraversalType() {
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
        logger.debug(String.format("TraverseReadsNano.traverse Covered dataset is %s", dataProvider));

        if( !dataProvider.hasReads() )
            throw new IllegalArgumentException("Unable to traverse reads; no read data is available.");

        nanoScheduler.setDebug(DEBUG);
        final TraverseReadsMap myMap = new TraverseReadsMap(walker);
        final TraverseReadsReduce myReduce = new TraverseReadsReduce(walker);

        T result = nanoScheduler.execute(aggregateMapData(dataProvider).iterator(), myMap, sum, myReduce);
        // TODO -- how do we print progress?
        //printProgress(dataProvider.getShard(), ???);

        return result;
    }

    private List<MapData> aggregateMapData(final ReadShardDataProvider dataProvider) {
        final ReadView reads = new ReadView(dataProvider);
        final ReadReferenceView reference = new ReadReferenceView(dataProvider);
        final ReadBasedReferenceOrderedView rodView = new ReadBasedReferenceOrderedView(dataProvider);

        final List<MapData> mapData = new ArrayList<MapData>();  // TODO -- need size of reads
        for ( final SAMRecord read : reads ) {
            final ReferenceContext refContext = ! read.getReadUnmappedFlag()
                    ? reference.getReferenceContext(read)
                    : null;

            // if the read is mapped, create a metadata tracker
            final RefMetaDataTracker tracker = read.getReferenceIndex() >= 0
                    ? rodView.getReferenceOrderedDataForRead(read)
                    : null;

            // update the number of reads we've seen
            dataProvider.getShard().getReadMetrics().incrementNumIterations();

            mapData.add(new MapData((GATKSAMRecord)read, refContext, tracker));
        }

        return mapData;
    }

    @Override
    public void printOnTraversalDone() {
        nanoScheduler.shutdown();
        super.printOnTraversalDone();
    }

    private class TraverseReadsReduce implements ReduceFunction<M, T> {
        final ReadWalker<M,T> walker;

        private TraverseReadsReduce(ReadWalker<M, T> walker) {
            this.walker = walker;
        }

        @Override
        public T apply(M one, T sum) {
            return walker.reduce(one, sum);
        }
    }

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

    private class TraverseReadsMap implements MapFunction<MapData, M> {
        final ReadWalker<M,T> walker;

        private TraverseReadsMap(ReadWalker<M, T> walker) {
            this.walker = walker;
        }

        @Override
        public M apply(final MapData data) {
            if ( ! walker.isDone() ) {
                final boolean keepMeP = walker.filter(data.refContext, data.read);
                if (keepMeP) {
                    return walker.map(data.refContext, data.read, data.tracker);
                }
            }

            return null; // TODO -- what should we return in the case where the walker is done or the read is filtered?
        }
    }
}
