package org.broadinstitute.sting.gatk.dataSources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;

/**
 *
 * User: aaron
 * Date: Apr 10, 2009
 * Time: 5:03:13 PM
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
 * @date Apr 10, 2009
 * <p/>
 * Class ReadShard
 * <p/>
 * A class for sharded reads.
 */
public class ReadShard implements Shard {

    // the count of the reads we want to copy off
    private int size = 0;

    // this is going to get gross
    private final ReadShardStrategy str;

    /**
     * create a read shard, given a read size
     *
     * @param size
     */
    public ReadShard(int size) {
        this.str = null;
        this.size = size;
    }

    /**
     * create a read shard, given a read size
     *
     * @param size
     */
    ReadShard(ReadShardStrategy caller, int size) {
        this.str = caller;
        this.size = size;
    }

    /** @return the genome location represented by this shard */
    public GenomeLoc getGenomeLoc() {
        throw new UnsupportedOperationException("Reads based sharding isn't genome loc aware");
    }

    /** @return the genome location represented by this shard */
    public int getSize() {
        return size;
    }


    /**
     * what kind of shard do we return
     *
     * @return ShardType, indicating the type
     */
    public ShardType getShardType() {
        return ShardType.READ;
    }
}
