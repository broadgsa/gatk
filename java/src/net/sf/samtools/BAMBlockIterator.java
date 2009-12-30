package net.sf.samtools;

import net.sf.samtools.util.BinaryCodec;
import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.StringLineReader;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SeekableStream;

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Arrays;
import java.nio.channels.FileChannel;
import java.io.IOException;
import java.io.FileInputStream;
import java.io.File;

/**
 * Walks over a BAM file, discovering and returning the starting location of each block
 * in chunk format.
 *
 * @author mhanna
 * @version 0.1
 */
public class BAMBlockIterator implements CloseableIterator<Block> {
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
     * @param file stream File to use when accessing the BAM.
     */
    public BAMBlockIterator(File file) throws IOException {
        FileInputStream inputStream = new FileInputStream(file);
        this.blockReader = new BlockReader(inputStream);
    }

    public void close() {
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
    public Block next() {
        if(!hasNext())
            throw new NoSuchElementException("No next chunk is available.");

        Block block = null;
        try {
            block = blockReader.getBlockAt(position);
            position = block.position + block.compressedBlockSize;
        }
        catch(IOException ex) {
            throw new SAMException("Unable to completely read chunk at end of file.", ex);            
        }
        return block;
    }

    /**
     * Remove a chunk from the iterator.
     * @throws UnsupportedOperationException always.
     */
    public void remove() {
        throw new UnsupportedOperationException("BAMChunkIterator: Cannot remove a BAM chunk.");
    }
}
