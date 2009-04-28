package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;

import java.util.Iterator;

/**
 *
 * User: aaron
 * Date: Apr 14, 2009
 * Time: 1:34:28 PM
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
 * A descriptions should go here. Blame aaron if it's missing.
 */
public class ReadShardStrategy implements ShardStrategy {

    // do we use unmapped reads in the sharding strategy
    private boolean unMappedReads = true;

    // our read bucket size, default
    protected long readCount = 100000L;

    // our sequence dictionary
    final private SAMSequenceDictionary dic;

    /**
     * the default constructor
     * @param dic the dictionary
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
        return true;
    }

    public Shard next() {
        return new ReadShard((int)readCount);  //To change body of implemented methods use File | Settings | File Templates.
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
}
