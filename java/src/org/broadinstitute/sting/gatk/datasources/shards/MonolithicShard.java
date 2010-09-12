package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GATKException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.ReadMetrics;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;

import java.util.List;

/**
 * A single, monolithic shard bridging all available data.
 * @author mhanna
 * @version 0.1
 */
public class MonolithicShard implements Shard {
    /**
     * Reads data, if applicable.
     */
    private final SAMDataSource readsDataSource;

    /**
     * What type of MonolithicShard is this?  Read or locus?
     */
    private final ShardType shardType;

    /**
     * Locations.  For the monolithic shard, should be a list of all available contigs in the reference.
     */
    private final List<GenomeLoc> locs;

    /**
     * Statistics about which reads in this shards were used and which were filtered away.
     */
    private final ReadMetrics readMetrics = new ReadMetrics();

    /**
     * Creates a new monolithic shard of the given type.
     * @param shardType Type of the shard.  Must be either read or locus; cannot be intervalic.
     * @param locs Intervals that this monolithic shard should process. 
     */
    public MonolithicShard(SAMDataSource readsDataSource, ShardType shardType, List<GenomeLoc> locs) {
        this.readsDataSource = readsDataSource;
        if(shardType != ShardType.LOCUS && shardType != ShardType.READ)
            throw new GATKException("Invalid shard type for monolithic shard: " + shardType);
        this.shardType = shardType;
        this.locs = locs;
    }

    /**
     * Closes the shard, tallying and incorporating read data.
     */
    @Override
    public void close() {
        readsDataSource.incorporateReadMetrics(readMetrics);
    }

    /**
     * Returns null, indicating that (in this case) the entire genome is covered.
     * @return null.
     */
    public List<GenomeLoc> getGenomeLocs() {
        return locs;
    }

    /**
     * Reports the type of monolithic shard.
     * @return Type of monolithic shard.
     */
    @Override
    public ShardType getShardType() {
        return shardType;
    }

    /**
     * Gets key read validation and filtering properties.
     * @return set of read properties associated with this shard.
     */
    @Override
    public ReadProperties getReadProperties() {
        return readsDataSource.getReadsInfo();
    }

    /**
     * Retrieves a storage space of metrics about number of reads included, filtered, etc.
     * @return Storage space for metrics.
     */
    @Override
    public ReadMetrics getReadMetrics() {
        return readMetrics;
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
