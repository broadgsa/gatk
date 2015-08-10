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
import htsjdk.samtools.SAMRecordCoordinateComparator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.datasources.providers.ReadShardDataProvider;
import org.broadinstitute.gatk.engine.datasources.providers.ReadView;
import org.broadinstitute.gatk.engine.datasources.reads.Shard;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.ReadPairWalker;
import org.broadinstitute.gatk.engine.walkers.Requires;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Traverse over a collection of read pairs, assuming that a given shard will contain all pairs.
 *
 * @author mhanna
 * @version 0.1
 */
@Requires({DataSource.REFERENCE})
public class TraverseReadPairs<M,T> extends TraversalEngine<M,T, ReadPairWalker<M,T>,ReadShardDataProvider> {

    /** our log, which we want to capture anything from this class */
    protected static final Logger logger = Logger.getLogger(TraverseReadPairs.class);

    @Override
    public String getTraversalUnits() {
        return "read pairs";
    }

    /**
     * Traverse by reads, given the data and the walker
     *
     * @param walker the walker to execute over
     * @param sum    of type T, the return from the walker
     *
     * @return the result type T, the product of all the reduce calls
     */
    public T traverse(ReadPairWalker<M, T> walker,
                      ReadShardDataProvider dataProvider,
                      T sum) {
        logger.debug(String.format("TraverseReadsPairs.traverse Covered dataset is %s", dataProvider));

        if( !dataProvider.hasReads() )
            throw new IllegalArgumentException("Unable to traverse reads; no read data is available.");

        ReadView reads = new ReadView(dataProvider);
        List<SAMRecord> pairs = new ArrayList<SAMRecord>();

        boolean done = walker.isDone();
        for(SAMRecord read: reads) {
            if ( done ) break;
            dataProvider.getShard().getReadMetrics().incrementNumReadsSeen();

            if(pairs.size() == 0 || pairs.get(0).getReadName().equals(read.getReadName())) {
                // If this read name is the same as the last, accumulate it.
                pairs.add(read);
            }
            else {
                // Otherwise, walk over the accumulated list, then start fresh with the new read.
                sum = walkOverPairs(walker,dataProvider.getShard(),pairs,sum);
                pairs.clear();
                pairs.add(read);

                printProgress(null);
            }

            done = walker.isDone();
        }

        // If any data was left in the queue, process it.
        if(pairs.size() > 0)
            sum = walkOverPairs(walker,dataProvider.getShard(),pairs,sum);

        return sum;
    }

    /**
     * Filter / map / reduce over a single pair.
     * @param walker The walker.
     * @param shard The shard currently being processed.
     * @param reads The reads in the pair.
     * @param sum The accumulator.
     * @return The accumulator after application of the given read pairing.
     */
    private T walkOverPairs(ReadPairWalker<M,T> walker, Shard shard, List<SAMRecord> reads, T sum) {
        // update the number of reads we've seen
        shard.getReadMetrics().incrementNumIterations();

        // Sort the reads present in coordinate order.
        Collections.sort(reads,new SAMRecordCoordinateComparator());

        final boolean keepMeP = walker.filter(reads);
        if (keepMeP) {
            M x = walker.map(reads);
            sum = walker.reduce(x, sum);
        }

        return sum;
    }
}
