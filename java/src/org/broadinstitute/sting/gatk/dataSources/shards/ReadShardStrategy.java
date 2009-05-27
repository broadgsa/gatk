package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;

import java.util.Iterator;

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
 * @date Apr 14, 2009
 * <p/>
 * Class ReadShardStrategy
 * <p/>
 * The sharding strategy for reads using a simple counting mechanism.  Each read shard
 * has a specific number of reads (default to 100K) which is configured in the constructor.
 */
public class ReadShardStrategy implements ShardStrategy {

    // do we use unmapped reads in the sharding strategy
    private boolean unMappedReads = true;

    // our read bucket size, default
    protected long readCount = 100000L;

    // our sequence dictionary
    final private SAMSequenceDictionary dic;

    // our hasnext flag
    boolean hasNext = true;

    /**
     * the default constructor
     * @param dic the sequence dictionary to use
     * @param size the read count to iterate over
     */
    ReadShardStrategy(SAMSequenceDictionary dic, long size) {
        this.dic = dic;
        readCount = size;    
    }

    /**
     * do we have another read shard?
     * @return
     */
    public boolean hasNext() {
        return hasNext;
    }

    public Shard next() {
        return new ReadShard((int)readCount, this);
    }

    public void remove() {
        throw new UnsupportedOperationException("Remove not supported");
    }

    public Iterator<Shard> iterator() {
        return this;
    }

    /**
     * set the next shards size
     *
     * @param size adjust the next size to this
     */
    public void adjustNextShardSize(long size) {
        readCount = size;
    }


    /**
     * this function is a work-around for the fact that
     * we don't know when we're out of reads until the SAM data source
     * tells us so.  
     */
    public void signalDone() {
        hasNext = false;    
    }
}
