package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;


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
 * <p/>
 * Class IntervalShard
 * <p/>
 * the base interval shard.  All interval shards are generally the same,
 * but must return their ShardType individually.  
 */
public class IntervalShard implements Shard {

    /** a collection of genomic locations to interate over */
    private GenomeLoc mSet;
    private Shard.ShardType mType = Shard.ShardType.LOCUS_INTERVAL;

    IntervalShard(GenomeLoc myLocation, Shard.ShardType intervalType) {
        if (intervalType != Shard.ShardType.LOCUS_INTERVAL && intervalType != Shard.ShardType.READ_INTERVAL)
            throw new IllegalArgumentException("The specified interval type must be either LOCUS_INTERVAL or READ_INTERVAL");
        mType = intervalType;
        mSet = myLocation.clone();
    }

    /** @return the genome location represented by this shard */
    public GenomeLoc getGenomeLoc() {
        return mSet;
    }

    /**
     * returns the type of shard, READ
     *
     * @return READ, indicating the shard type
     */
    public Shard.ShardType getShardType() {
        return mType;
    }

    /**
     * String representation of this shard.
     * @return A string representation of the boundaries of this shard.
     */
    @Override
    public String toString() {
        return mSet.toString();
    }    
}
