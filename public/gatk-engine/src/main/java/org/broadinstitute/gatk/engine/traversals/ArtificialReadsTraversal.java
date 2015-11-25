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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.datasources.providers.ShardDataProvider;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.sam.ArtificialPatternedSAMIterator;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;


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
public class ArtificialReadsTraversal<M,T> extends TraversalEngine<M,T,Walker<M,T>,ShardDataProvider> {

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

    @Override
    public String getTraversalUnits() {
        return "reads";
    }

    /**
     * Traverse by reads, given the data and the walker
     *
     * @param walker       the walker to traverse with
     * @param dataProvider the provider of the reads data
     * @param sum          the value of type T, specified by the walker, to feed to the walkers reduce function
     *
     * @return the reduce variable of the read walker
     */
    public T traverse( Walker<M, T> walker,
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

            // an array of characters that represent the reference
            ReferenceContext refSeq = null;

            final boolean keepMeP = readWalker.filter(refSeq, (GATKSAMRecord) read);
            if (keepMeP) {
                M x = readWalker.map(refSeq, (GATKSAMRecord) read, null);  // TODO: fix me at some point, it would be nice to fake out ROD data too
                sum = readWalker.reduce(x, sum);
            }
        }
        return sum;
    }
}
