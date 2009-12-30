package net.sf.samtools;

import net.sf.samtools.util.BlockCompressedStreamConstants;

import java.io.IOException;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;

/**
 * Read an individual block from the BGZF file.
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
     * @param inputStream InputStream from which to read.
     */
    public BlockReader(final FileInputStream inputStream) {
        this.channel = inputStream.getChannel();
        this.buffer = ByteBuffer.allocateDirect(BlockCompressedStreamConstants.MAX_COMPRESSED_BLOCK_SIZE);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
    }

    /**
     * Closes the block reader's channel.
     * @throws IOException On failure to close channel.
     */
    public void close() throws IOException {
        this.channel.close();
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
    public Block getBlockAt(long position) throws IOException {
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
        final int compressedBlockSize = (buffer.getShort() & 0xffff) + 1;
        if (compressedBlockSize < BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH || compressedBlockSize > buffer.capacity()) {
            throw new IOException("Unexpected compressed block length: " + compressedBlockSize);
        }

        // Read the uncompressed block size
        buffer.rewind();
        buffer.limit(4);
        channel.read(buffer,position+compressedBlockSize-4);
        buffer.flip();
        final int uncompressedBlockSize = buffer.getInt();

        return new Block(position,compressedBlockSize,uncompressedBlockSize);
    }
}
