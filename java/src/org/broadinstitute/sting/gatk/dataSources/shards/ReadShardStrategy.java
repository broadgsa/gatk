package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;

import java.util.Iterator;

/**
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
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
