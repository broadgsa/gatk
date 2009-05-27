package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;

import java.util.List;

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
 * @version 1.0
 * @date Apr 6, 2009
 * <p/>
 * Class AdaptiveShard
 * <p/>
 * allows you to change the sharding length as you traverse
 */
class LinearLocusShardStrategy extends LocusShardStrategy {

    // default the next size to 100,000
    private long nextShardSize = 100000;

    /**
     * the constructor, taking a seq dictionary to parse out contigs
     *
     * @param dic the seq dictionary
     */
    LinearLocusShardStrategy(SAMSequenceDictionary dic, long startSize) {
        super(dic);
        this.nextShardSize = startSize;
    }

    /**
     * the constructor, constructing from another Locus strategy
     *
     * @param strat the shatter to convert from
     */
    LinearLocusShardStrategy(LocusShardStrategy strat) {
        super(strat);
        this.nextShardSize = strat.nextShardSize();
    }


    /**
     * The constructor, for a genomic list, start size, and a reference dictionary 
     * @param dic the reference dictionary
     * @param startSize the starting size of the shard
     * @param lst locations to iterate from
     */
    LinearLocusShardStrategy(SAMSequenceDictionary dic, long startSize, GenomeLocSortedSet lst) {
        super(dic, lst);
        this.nextShardSize = startSize;
    }


    /**
     * set the next shards size
     *
     * @param size adjust the next size to this
     */
    public void adjustNextShardSize(long size) {
        nextShardSize = size;
    }

    /**
     * This is how the various shards strategies implements their approach
     *
     * @return the next shard size
     */
    protected long nextShardSize() {
        return nextShardSize;
    }

}
