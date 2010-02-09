package org.broadinstitute.sting.gatk.datasources.shards;

import net.sf.samtools.Chunk;
import net.sf.samtools.SAMFileReader2;

import java.util.List;
import java.util.Map;
import java.util.Collections;

/**
 * Expresses a shard of read data in block format.
 *
 * @author mhanna
 * @version 0.1
 */
public class BlockDelimitedReadShard extends ReadShard implements BAMFormatAwareShard {
    /**
     * The reader associated with this shard.
     */
    private final SAMFileReader2 reader;
    
    /**
     * The list of chunks to retrieve when loading this shard.
     */
    private final List<Chunk> chunks;

    public BlockDelimitedReadShard(SAMFileReader2 reader, List<Chunk> chunks) {
        this.reader = reader;
        this.chunks = chunks;
    }

    /**
     * Get the list of chunks delimiting this shard.
     * @return a list of chunks that contain data for this shard.
     */
    public Map<SAMFileReader2,List<Chunk>> getChunks() {
        return Collections.singletonMap(reader,chunks);
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
