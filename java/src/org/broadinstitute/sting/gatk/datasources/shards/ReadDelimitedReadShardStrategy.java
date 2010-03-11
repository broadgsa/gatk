package org.broadinstitute.sting.gatk.datasources.shards;

import java.util.Iterator;

/**
 * A shard strategy that breaks up shards based on how many reads are
 * in each.
 *
 * @author mhanna
 * @version 0.1
 */
public class ReadDelimitedReadShardStrategy extends ReadShardStrategy {
    // our read bucket size, default
    protected long readCount = 1000L;

    // our hasnext flag
    boolean hasNext = true;

    // our limiting factor
    long limitedSize = -1;
    boolean stopDueToLimitingFactor = false;

    /**
     * the default constructor
     * @param size the read count to iterate over
     * @param limitedSize limit the shard to this length
     */
    ReadDelimitedReadShardStrategy(long size, long limitedSize) {
        readCount = size;
        this.limitedSize = limitedSize;
    }

    /**
     * do we have another read shard?
     * @return
     */
    public boolean hasNext() {
        if (stopDueToLimitingFactor) {
            return false;
        }
        return hasNext;
    }

    public Shard next() {
        if (limitedSize > 0) {
            if (limitedSize > readCount) {
                limitedSize = limitedSize - readCount;
            }
            else {
                readCount = limitedSize;
                limitedSize = 0;
                stopDueToLimitingFactor = true;
            }
        }

        return new ReadDelimitedReadShard(this,(int)readCount);
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
