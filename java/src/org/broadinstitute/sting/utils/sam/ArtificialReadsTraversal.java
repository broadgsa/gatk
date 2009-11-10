package org.broadinstitute.sting.utils.sam;

import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.gatk.traversals.TraversalStatistics;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.apache.log4j.Logger;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;


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

/**
 * @author aaron
 *
 * this class acts as a fake reads traversal engine for testing out reads based traversals.
 */
public class ArtificialReadsTraversal extends TraversalEngine {

    public int startingChr = 1;
    public int endingChr = 5;
    public int readsPerChr = 100;
    public int unMappedReads = 1000;
    private int DEFAULT_READ_LENGTH = ArtificialSAMUtils.DEFAULT_READ_LENGTH;
    private ArtificialPatternedSAMIterator iter;
    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(ArtificialReadsTraversal.class);

    /** Creates a new, uninitialized ArtificialReadsTraversal */
    public ArtificialReadsTraversal() {
    }

    // what read ordering are we using
    private ArtificialPatternedSAMIterator.PATTERN readOrder = ArtificialPatternedSAMIterator.PATTERN.IN_ORDER_READS;


    /**
     * set the read ordering of the reads given to the walker
     *
     * @param readOrdering
     */
    public void setReadOrder( ArtificialPatternedSAMIterator.PATTERN readOrdering ) {
        readOrder = readOrdering;
    }

    /**
     * Traverse by reads, given the data and the walker
     *
     * @param walker       the walker to traverse with
     * @param shard        the shard, specifying the range of data to iterate over
     * @param dataProvider the provider of the reads data
     * @param sum          the value of type T, specified by the walker, to feed to the walkers reduce function
     * @param <M>          the map type of the walker
     * @param <T>          the reduce type of the walker
     *
     * @return the reduce variable of the read walker
     */
    public <M, T> T traverse( Walker<M, T> walker,
                              Shard shard,
                              ShardDataProvider dataProvider,
                              T sum ) {

        if (!( walker instanceof ReadWalker ))
            throw new IllegalArgumentException("Walker isn't a read walker!");

        ReadWalker<M, T> readWalker = (ReadWalker<M, T>) walker;
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(( endingChr - startingChr ) + 1, startingChr, readsPerChr + DEFAULT_READ_LENGTH);
        iter = new ArtificialPatternedSAMIterator(this.startingChr,
                this.endingChr,
                this.readsPerChr,
                this.unMappedReads,
                header,
                this.readOrder);

        // while we still have more reads
        for (SAMRecord read : iter) {

            // our alignment context
            AlignmentContext alignment = null;

            // an array of characters that represent the reference
            char[] refSeq = null;

            // update the number of reads we've seen
            TraversalStatistics.nRecords++;

            final boolean keepMeP = readWalker.filter(refSeq, read);
            if (keepMeP) {
                M x = readWalker.map(refSeq, read);
                sum = readWalker.reduce(x, sum);
            }

            if (alignment != null) { printProgress(TRAVERSAL_TYPE.READ, alignment.getLocation()); }
        }
        return sum;
    }

    /**
     * Temporary override of printOnTraversalDone.
     * TODO: Add some sort of TE.getName() function once all TraversalEngines are ported.
     *
     * @param sum Result of the computation.
     * @param <T> Type of the result.
     */
    public <T> void printOnTraversalDone( T sum ) {
        printOnTraversalDone(TRAVERSAL_TYPE.READ, sum);
    }


}
