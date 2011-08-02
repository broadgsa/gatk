package org.broadinstitute.sting.gatk.datasources.reads;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Create a giant shard representing all the data in the input BAM(s).
 *
 * @author mhanna
 * @version 0.1
 */
public class MonolithicShardStrategy implements ShardStrategy {
    /**
     * The single shard associated with this sharding strategy.
     */
    private MonolithicShard shard;

    /**
     * Create a new shard strategy for shards of the given type.
     * @param shardType The shard type.
     */
    public MonolithicShardStrategy(final GenomeLocParser parser, final SAMDataSource readsDataSource, final Shard.ShardType shardType, final List<GenomeLoc> region) {
        shard = new MonolithicShard(parser,readsDataSource,shardType,region);
    }

    /**
     * Convenience for using in a foreach loop.  Will NOT create a new, reset instance of the iterator;
     * will only return another copy of the active iterator. 
     * @return A copy of this.
     */
    public Iterator<Shard> iterator() {
        return this;
    }

    /**
     * Returns true if the monolithic shard has not yet been consumed, or false otherwise.
     * @return True if shard has been consumed, false otherwise.
     */
    public boolean hasNext() {
        return shard != null;
    }

    /**
     * Returns the monolithic shard if it has not already been retrieved.
     * @return The monolithic shard.
     * @throws NoSuchElementException if no such data exists.
     */
    public Shard next() {
        if(shard == null)
            throw new NoSuchElementException("Monolithic shard has already been retrived.");

        Shard working = shard;
        shard = null;
        return working;
    }

    /**
     * Mandated by the interface, but is unsupported in this context.  Will throw an exception always.
     */
    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a shard strategy");
    }

    /**
     * Mandated by the interface, but is unsupported in this context.  Will throw an exception always.
     * @param size adjust the next size to this
     */
    public void adjustNextShardSize( long size ) {
        throw new UnsupportedOperationException("Cannot adjust the next size of a monolithic shard; there will be no next shard.");
    }

}

