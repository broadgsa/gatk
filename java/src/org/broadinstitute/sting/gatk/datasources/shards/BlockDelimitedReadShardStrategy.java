package org.broadinstitute.sting.gatk.datasources.shards;

import net.sf.samtools.Chunk;
import net.sf.samtools.BAMFileHeaderLoader;
import net.sf.samtools.BAMChunkIterator;
import net.sf.samtools.BAMBlockIterator;

import java.util.*;
import java.io.File;
import java.io.IOException;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;

/**
 * A read shard strategy that delimits based on the number of
 * blocks in the BAM file.
 *
 * @author mhanna
 * @version 0.1
 */
public class BlockDelimitedReadShardStrategy extends ReadShardStrategy {
    /**
     * Number of blocks in a given shard.
     */
    protected int blockCount = 100;

    /**
     * The actual chunks streaming into the file.
     */
    private final BAMChunkIterator chunkIterator;

    /**
     * The data backing the next chunks to deliver to the traversal engine.
     */
    private final List<Chunk> nextChunks;

    /**
     * Create a new read shard strategy, loading read shards from the given BAM file.
     * @param dataSource Data source from which to load shards.
     */
    public BlockDelimitedReadShardStrategy(SAMDataSource dataSource) {
        if(dataSource.getReadsInfo().getReadsFiles().size() > 1)
            throw new UnsupportedOperationException("Experimental sharding only works with a single BAM at the moment.");
        File bamFile = dataSource.getReadsInfo().getReadsFiles().get(0);
        try {
            Chunk headerLocation = BAMFileHeaderLoader.load(bamFile).getLocation();
            chunkIterator = new BAMChunkIterator(new BAMBlockIterator(bamFile),Arrays.asList(BAMFileHeaderLoader.preambleLocation,headerLocation));
        }
        catch(IOException ex) {
            throw new StingException("Unable to open BAM file for sharding.");
        }
        nextChunks = new ArrayList<Chunk>();

        advance();
    }

    /**
     * do we have another read shard?
     * @return True if any more data is available.  False otherwise.
     */
    public boolean hasNext() {
        return nextChunks.size() > 0;
    }

    /**
     * Retrieves the next shard, if available.
     * @return The next shard, if available.
     * @throws NoSuchElementException if no such shard is available.
     */
    public Shard next() {
        if(!hasNext())
            throw new NoSuchElementException("No such element available: SAM reader has arrived at last shard.");
        Shard shard = new BlockDelimitedReadShard(Collections.unmodifiableList(new ArrayList<Chunk>(nextChunks)));
        advance();
        return shard;
    }

    /**
     * @throws UnsupportedOperationException always.
     */
    public void remove() {
        throw new UnsupportedOperationException("Remove not supported");
    }

    /**
     * Convenience method for using ShardStrategy in an foreach loop.
     * @return A iterator over shards.
     */
    public Iterator<Shard> iterator() {
        return this;
    }

    private void advance() {
        nextChunks.clear();
        int chunksCopied = 0;

        while(chunksCopied++ < blockCount && chunkIterator.hasNext())
            nextChunks.add(chunkIterator.next());
    }
}
