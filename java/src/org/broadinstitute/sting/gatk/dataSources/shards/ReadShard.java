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
 *         <p/>
 *         ReadShard
 *         <p/>
 *         the base class for read shards.
 */
public class ReadShard implements Shard {

    // the count of the reads we want to copy off
    private int size = 0;

    /**
     * our tie in for the shard strategy.  This allows us to signal to the shard
     * strategy that we've finished process, so it can indicate that we're out of reads
     */
    private final ReadShardStrategy str;

    // the reference back to our read shard strategy
    private final ReadShardStrategy strat;

    /**
     * create a read shard, given a read size
     *
     * @param size
     */
    ReadShard(int size, ReadShardStrategy strat) {
        this.str = null;
        this.size = size;
        this.strat = strat;
    }

    /**
     * create a read shard, given a read size
     *
     * @param size
     */
    ReadShard(ReadShardStrategy caller, int size, ReadShardStrategy strat) {
        this.str = caller;
        this.size = size;
        this.strat = strat;
    }

    /** @return the genome location represented by this shard */
    public GenomeLoc getGenomeLoc() {
        throw new UnsupportedOperationException("ReadShard isn't genome loc aware");
    }

    /** @return the genome location represented by this shard */
    public int getSize() {
        return size;
    }

    /**
     * this method is used as a backend, to signal to the sharding strategy that we've
     * finished processing.  When we move to a more read-aware bam system this method could disappear. 
     */
    public void signalDone() {
        strat.signalDone();
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
