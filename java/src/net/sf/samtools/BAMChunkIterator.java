package net.sf.samtools;

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.nio.channels.FileChannel;
import java.io.IOException;

/**
 * Walks over a BAM file, discovering and returning the starting location of each block
 * in chunk format.
 *
 * @author mhanna
 * @version 0.1
 */
public class BAMChunkIterator implements Iterator<Chunk> {
    /**
     * File channel from which to read chunks.  
     */
    private BlockReader blockReader;

    /**
     * What is the current position of this block within the BAM file?
     */
    private long position = 0;

    /**
     * Iterate through the BAM chunks in a file.
     * @param channel File channel to use when accessing the BAM.
     */
    public BAMChunkIterator(FileChannel channel) {
        this.blockReader = new BlockReader(channel);
    }

    /**
     * Are there other chunks to retrieve in this file?
     * @return True if other chunks are available, false otherwise.
     */
    public boolean hasNext() {
        try {
            return !blockReader.eof(position);
        } catch(IOException ex) {
            throw new SAMException("Unable to check file for a next BAM record", ex);
        }
    }

    /**
     * Get the next chunk from the iterator.
     * @return The next chunk.
     * @throw NoSuchElementException if no next chunk is available.
     */
    public Chunk next() {
        if(!hasNext())
            throw new NoSuchElementException("No next chunk is available.");

        Chunk chunk = null;
        try {
            chunk = blockReader.getChunkAt(position);
            position = (chunk.getChunkEnd() >> 16) + 1;
        }
        catch(IOException ex) {
            throw new SAMException("Unable to completely read chunk at end of file.", ex);            
        }
        return chunk;
    }

    /**
     * Remove a chunk from the iterator.
     * @throws UnsupportedOperationException always.
     */
    public void remove() {
        throw new UnsupportedOperationException("BAMChunkIterator: Cannot remove a BAM chunk.");
    }
}
