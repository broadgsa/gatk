package org.broadinstitute.sting.gatk.datasources.shards;

import net.sf.samtools.Chunk;

import java.util.List;

/**
 * A common interface for shards that natively understand the BAM format.
 *
 * @author mhanna
 * @version 0.1
 */
public interface BAMFormatAwareShard extends Shard {
    /**
     * Get the list of chunks delimiting this shard.
     * @return a list of chunks that contain data for this shard.
     */
    public List<Chunk> getChunks();    
}
