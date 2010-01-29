package org.broadinstitute.sting.gatk.datasources.shards;

import net.sf.samtools.Chunk;

import java.util.List;

/**
 * Expresses a shard of read data in block format.
 *
 * @author mhanna
 * @version 0.1
 */
public class BlockDelimitedReadShard extends ReadShard implements BAMFormatAwareShard {
    /**
     * The list of chunks to retrieve when loading this shard.
     */
    private final List<Chunk> chunks;

    public BlockDelimitedReadShard(List<Chunk> chunks) {
        this.chunks = chunks;
    }

    /**
     * Get the list of chunks delimiting this shard.
     * @return a list of chunks that contain data for this shard.
     */
    public List<Chunk> getChunks() {
        return chunks;
    }

    /**
     * String representation of this shard.
     * @return A string representation of the boundaries of this shard.
     */    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for(Chunk chunk : chunks) {
            sb.append(chunk);
            sb.append(' ');
        }
        return sb.toString();
    }
}
