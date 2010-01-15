package org.broadinstitute.sting.gatk.datasources.shards;

/**
 * A read shard delimited by an actual read count, rather than blocks or any other
 * physical mapping of the BAM file.
 *
 * @author mhanna
 * @version 0.1
 */
public class ReadDelimitedReadShard extends ReadShard {
    // the count of the reads we want to copy off
    private int size = 0;

    /**
     * our tie in for the shard strategy.  This allows us to signal to the shard
     * strategy that we've finished process, so it can indicate that we're out of reads
     */
    private final ReadDelimitedReadShardStrategy strat;

    /**
     * create a read shard, given a read size
     * @param strat The sharding strategy used to create this shard.
     * @param size Size of the shard, in reads.
     */
    ReadDelimitedReadShard(ReadDelimitedReadShardStrategy strat, int size) {
        this.size = size;
        this.strat = strat;
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
     * String representation of this shard.
     * @return A string representation of the boundaries of this shard.
     */
    @Override
    public String toString() {
        return String.format("%d reads", size);
    }
}
