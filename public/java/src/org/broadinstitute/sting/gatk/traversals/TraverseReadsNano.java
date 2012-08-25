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
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.datasources.reads.ReadShard;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.nanoScheduler.MapFunction;
import org.broadinstitute.sting.utils.nanoScheduler.NanoScheduler;
import org.broadinstitute.sting.utils.nanoScheduler.ReduceFunction;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

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
    final NanoScheduler<SAMRecord, M, T> nanoScheduler;

    public TraverseReadsNano(int nThreads) {
        final int bufferSize = ReadShard.getReadBufferSize() + 1; // actually has 1 more than max
        final int mapGroupSize = bufferSize / 10 + 1;
        nanoScheduler = new NanoScheduler<SAMRecord, M, T>(bufferSize, mapGroupSize, nThreads);
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

        if ( dataProvider.hasReferenceOrderedData() )
            throw new ReviewedStingException("Parallel read walkers currently don't support access to reference ordered data");

        final ReadView reads = new ReadView(dataProvider);
        final ReadReferenceView reference = new NotImplementedReadReferenceView(dataProvider);
        final ReadBasedReferenceOrderedView rodView = new ReadBasedReferenceOrderedView(dataProvider);

        nanoScheduler.setDebug(DEBUG);
        final TraverseReadsMap myMap = new TraverseReadsMap(reads, reference, rodView, walker);
        final TraverseReadsReduce myReduce = new TraverseReadsReduce(walker);

        T result = nanoScheduler.execute(reads.iterator().iterator(), myMap, sum, myReduce);
        // TODO -- how do we print progress?
        //printProgress(dataProvider.getShard(), ???);

        return result;
    }

    @Override
    public void printOnTraversalDone() {
        nanoScheduler.shutdown();
        super.printOnTraversalDone();    //To change body of overridden methods use File | Settings | File Templates.
    }

    private static class NotImplementedReadReferenceView extends ReadReferenceView {
        private NotImplementedReadReferenceView(ShardDataProvider provider) {
            super(provider);
        }

        @Override
        protected byte[] getReferenceBases(SAMRecord read) {
            throw new ReviewedStingException("Parallel read walkers don't support accessing reference yet");
        }

        @Override
        protected byte[] getReferenceBases(GenomeLoc genomeLoc) {
            throw new ReviewedStingException("Parallel read walkers don't support accessing reference yet");
        }
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

    private class TraverseReadsMap implements MapFunction<SAMRecord, M> {
        final ReadView reads;
        final ReadReferenceView reference;
        final ReadBasedReferenceOrderedView rodView;
        final ReadWalker<M,T> walker;

        private TraverseReadsMap(ReadView reads, ReadReferenceView reference, ReadBasedReferenceOrderedView rodView, ReadWalker<M, T> walker) {
            this.reads = reads;
            this.reference = reference;
            this.rodView = rodView;
            this.walker = walker;
        }

        @Override
        public M apply(final SAMRecord read) {
            if ( ! walker.isDone() ) {
                // ReferenceContext -- the reference bases covered by the read
                final ReferenceContext refContext = ! read.getReadUnmappedFlag() && reference != null
                        ? reference.getReferenceContext(read)
                        : null;

                // update the number of reads we've seen
                //dataProvider.getShard().getReadMetrics().incrementNumIterations();

                // if the read is mapped, create a metadata tracker
                final ReadMetaDataTracker tracker = read.getReferenceIndex() >= 0 ? rodView.getReferenceOrderedDataForRead(read) : null;

                final boolean keepMeP = walker.filter(refContext, (GATKSAMRecord) read);
                if (keepMeP) {
                    return walker.map(refContext, (GATKSAMRecord) read, tracker);
                }
            }

            return null; // TODO -- what should we return in the case where the walker is done or the read is filtered?
        }
    }
}
