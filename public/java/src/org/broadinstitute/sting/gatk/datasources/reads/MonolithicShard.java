package org.broadinstitute.sting.gatk.datasources.reads;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.List;

/**
 * A single, monolithic shard bridging all available data.
 * @author mhanna
 * @version 0.1
 */
public class MonolithicShard extends Shard {
    /**
     * Creates a new monolithic shard of the given type.
     * @param shardType Type of the shard.  Must be either read or locus; cannot be intervalic.
     * @param locs Intervals that this monolithic shard should process. 
     */
    public MonolithicShard(GenomeLocParser parser, SAMDataSource readsDataSource, ShardType shardType, List<GenomeLoc> locs) {
        super(parser, shardType, locs, readsDataSource, null, false);
        if(shardType != ShardType.LOCUS && shardType != ShardType.READ)
            throw new ReviewedStingException("Invalid shard type for monolithic shard: " + shardType);
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
