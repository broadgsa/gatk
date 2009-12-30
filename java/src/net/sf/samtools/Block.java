package net.sf.samtools;

/**
 * Represents the position of a block on disk.
 * 
 * @author mhanna
 * @version 0.1
 */
public class Block {
    public final long position;
    public final int compressedBlockSize;
    public final long uncompressedBlockSize;

    /**
     * Create a block, loading no data into memory.
     * @param position Position of this block on disk.s
     * @param compressedBlockSize Size of the block on disk; if compressedData is present, should match compressedData.length.
     * @param uncompressedBlockSize Size of the data in the block.
     */
    public Block(final long position, final int compressedBlockSize, final long uncompressedBlockSize) {
        this.position = position;
        this.compressedBlockSize = compressedBlockSize;
        this.uncompressedBlockSize = uncompressedBlockSize;
    }

    /**
     * Build a string representation of the block.
     * @return A string indicating position and size.
     */
    @Override
    public String toString() {
        return String.format("Block: pos = %d, compressed size = %d, uncompressed size = %d",position,compressedBlockSize,uncompressedBlockSize);
    }
}
