package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * A single, monolithic shard bridging all available data.
 * @author mhanna
 * @version 0.1
 */
public class MonolithicShard implements Shard {
    /**
     * What type of MonolithicShard is this?  Read or locus?
     */
    private ShardType shardType;

    /**
     * Creates a new monolithic shard of the given type.
     * @param shardType Type of the shard.  Must be either read or locus; cannot be intervalic.
     */
    public MonolithicShard(ShardType shardType) {
        if(shardType != ShardType.LOCUS && shardType != ShardType.READ)
            throw new StingException("Invalid shard type for monolithic shard: " + shardType);
        this.shardType = shardType;
    }

    /**
     * Returns null, indicating that (in this case) the entire genome is covered.
     * @return null.
     */
    public GenomeLoc getGenomeLoc() {
        return null;
    }

    /**
     * Reports the type of monolithic shard.
     * @return Type of monolithic shard.
     */
    public ShardType getShardType() {
        return shardType;
    }

    /**
     * String representation of this shard.
     * @return "entire genome".
     */
    @Override
    public String toString() {
        return "entire genome";    
    }
}
