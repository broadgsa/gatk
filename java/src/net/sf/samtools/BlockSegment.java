package net.sf.samtools;

import java.util.List;
import java.util.Collections;
import java.util.ArrayList;

/**
     * Represents a segment of a given block.  A block segment can easily be
 * converted or compared to either a block or a chunk.
 */
class BlockSegment {
    /**
     * Offset to the start of this block within the file.
     */
    private final long position;

    /**
     * Starting coordinate within this segment.
     */
    private final long blockStart;

    /**
     * Last coordinate within this segment.
     */
    private final long blockStop;

    /**
     * Create a new block segment from the given block.
     * @param block the block to use as the starting point for this segment.
     */
    public BlockSegment(Block block) {
        this(block.position,0,block.uncompressedBlockSize-1);
    }

    BlockSegment(final long position,final long blockStart,final long blockStop) {
        this.position = position;
        this.blockStart = blockStart;
        this.blockStop = blockStop;
    }

    /**
     * The size of the block segment.
     * @return Number of bytes represented by this block segment.
     */
    public long size() {
        return blockStop - blockStart + 1;
    }

    /**
     * Convert the given segment to a chunk.
     * @return the chunk equivalent of this block.
     */
    public Chunk toChunk() {
        return new Chunk(position << 16 & blockStart,position << 16 & blockStop);
    }

    /**
     * Does this block segment overlap the given chunk?
     * @param chunk The chunk to test.
     * @return true if the chunk is w/i the block; false otherwise.
     */
    public boolean overlaps(Chunk chunk) {
        Chunk segmentAsChunk = toChunk();

        // Chunk is completely encompasses the given block.
        if(chunk.getChunkStart() < segmentAsChunk.getChunkStart() &&
           chunk.getChunkEnd() > segmentAsChunk.getChunkEnd())
            return true;

        // Chunk starts w/i block.
        if(chunk.getChunkStart() >= segmentAsChunk.getChunkStart() &&
           chunk.getChunkStart() <= segmentAsChunk.getChunkEnd())
            return true;

        // Chunk ends w/i block.
        if(chunk.getChunkEnd() >= segmentAsChunk.getChunkStart() &&
           chunk.getChunkEnd() <= segmentAsChunk.getChunkEnd())
            return true;

        return false;
    }

    /**
     * Subtracts the given chunk from the block segment.
     * @param chunk The chunk to subtract.
     * @return The block segment(s) derived from subtracting the chunk.
     */
    public List<BlockSegment> minus(Chunk chunk) {
        // If no overlap exists, just return a new instance of the current block segment.
        if(!overlaps(chunk))
            return Collections.singletonList(this);

        List<BlockSegment> differenceSegments = new ArrayList<BlockSegment>();
        Chunk segmentAsChunk = toChunk();

        // If the segment begins before the chunk, build out a segment from the start of the block to just before
        // the start of the chunk or the end of the block segment, whichever comes first.
        if(segmentAsChunk.getChunkStart() < chunk.getChunkStart()) {
            long segmentStop = (position == chunk.getChunkStart()>>16) ? (chunk.getChunkStart()&0xFFFF)-1 : blockStop;
            differenceSegments.add(new BlockSegment(position,blockStart,segmentStop));
        }

        // If the segment ends after the chunk, build out a segment from either after the end of the chunk or the
        // start of the block segment, whichever comes last.
        if(segmentAsChunk.getChunkEnd() > chunk.getChunkEnd()) {
            long segmentStart = (position == chunk.getChunkEnd()>>16) ? (chunk.getChunkEnd()&0xFFFF)+1 : blockStart;
            differenceSegments.add(new BlockSegment(position,segmentStart,blockStop));
        }

        return differenceSegments;
    }
}
