package net.sf.samtools;

import net.sf.samtools.util.BlockCompressedStreamConstants;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;

/**
 * Read an individual block or chunk from the BGZF file.
 *
 * @author mhanna
 * @version 0.1
 */
public class BlockReader {
    /**
     * File channel from which to read.
     */
    private final FileChannel channel;

    /**
     * Store buffered information from the file temporarily.
     */
    private final ByteBuffer buffer;

    /**
     * Create a new block reader.  Block readers can operate independently on the same input file.
     * @param channel File channel from which to read.
     */
    public BlockReader(final FileChannel channel) {
        this.channel = channel;
        this.buffer = ByteBuffer.allocateDirect(BlockCompressedStreamConstants.MAX_COMPRESSED_BLOCK_SIZE);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
    }

    /**
     * Determine whether the given offset is at the end of the file.
     * @param position Position at which to start reading.  Must be at the beginning of a block.
     * @return True if we're at the end of the file (or at the end of the data stream), false otherwise.
     * @throws IOException If whether the position is past the end of the file can't be determined.
     */
    public boolean eof(long position) throws IOException {
        return (channel.size() - position) <= 0 ||
               (channel.size() - position) == BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length;
    }

    /**
     * Load the chunk information at a given position.
     * @param position Position at which the chunk begins.
     * @return Chunk of the BGZF file, starting at position.
     * @throws IOException if the chunk could not be read.
     */
    public Chunk getChunkAt(long position) throws IOException {
        buffer.rewind();
        buffer.limit(BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH);
        int count = channel.read(buffer,position);
        if (count == 0) {
            throw new IOException("Tried to read chunk past end of file");
        }
        if (count != BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH) {
            throw new IOException("Premature end of file");
        }
        buffer.flip();

        // Validate GZIP header
        if (buffer.get() != BlockCompressedStreamConstants.GZIP_ID1 ||
            buffer.get() != (byte)BlockCompressedStreamConstants.GZIP_ID2 ||
            buffer.get() != BlockCompressedStreamConstants.GZIP_CM_DEFLATE ||
            buffer.get() != BlockCompressedStreamConstants.GZIP_FLG) {
            throw new SAMFormatException("Invalid GZIP header");
        }
        // Skip MTIME, XFL, OS fields
        buffer.position(buffer.position() + 6);
        if (buffer.getShort() != BlockCompressedStreamConstants.GZIP_XLEN) {
            throw new SAMFormatException("Invalid GZIP header");
        }
        // Skip blocksize subfield intro
        buffer.position(buffer.position() + 4);
        // Read ushort
        final int totalBlockSize = (buffer.getShort() & 0xffff) + 1;
        if (totalBlockSize < BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH || totalBlockSize > buffer.capacity()) {
            throw new IOException("Unexpected compressed block length: " + totalBlockSize);
        }

        final long chunkStart = position << 16;
        final long chunkEnd = (position+totalBlockSize-1) << 16;

        return new Chunk(chunkStart,chunkEnd);
    }
}
