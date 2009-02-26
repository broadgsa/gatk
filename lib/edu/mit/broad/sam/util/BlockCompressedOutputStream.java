/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam.util;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.CRC32;
import java.util.zip.Deflater;

/**
 * Writer for a file that is a series of gzip blocks.  The caller just treats it as an
 * OutputStream, and under the covers a gzip block is written when the amount of uncompressed as-yet-unwritten
 * bytes reaches a threshold.  Note that the flush() method should not be called by client
 * unless you know what you're doing, because it forces a gzip block to be written even if the
 * number of buffered bytes has not reached threshold.  close(), on the other hand, must be called
 * when done writing in order to force the last gzip block to be written.
 */
public class BlockCompressedOutputStream
        extends OutputStream
{
    private final BinaryCodec codec;
    private final byte[] uncompressedBuffer = new byte[BlockCompressedStreamConstants.DEFAULT_UNCOMPRESSED_BLOCK_SIZE];
    private int numUncompressedBytes = 0;
    private final byte[] compressedBuffer =
            new byte[BlockCompressedStreamConstants.MAX_COMPRESSED_BLOCK_SIZE -
                    BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH];
    private final Deflater deflater = new Deflater(BlockCompressedStreamConstants.GZIP_CM_DEFLATE, true);
    private final CRC32 crc32 = new CRC32();
    private final byte[] singleByteArray = new byte[1];

    private int numberOfThrottleBacks = 0;

    public BlockCompressedOutputStream(final String filename) {
        codec = new BinaryCodec(filename, true);
    }

    public BlockCompressedOutputStream(final File file) {
        codec = new BinaryCodec(file, true);
    }

    @Override
    public void write(final byte[] bytes) throws IOException {
        write(bytes, 0, bytes.length);
    }

    @Override
    public void write(final byte[] bytes, int startIndex, int numBytes) throws IOException {
        assert(numUncompressedBytes < uncompressedBuffer.length);
        while (numBytes > 0) {
            final int bytesToWrite = Math.min(uncompressedBuffer.length - numUncompressedBytes, numBytes);
            System.arraycopy(bytes, startIndex, uncompressedBuffer, numUncompressedBytes, bytesToWrite);
            numUncompressedBytes += bytesToWrite;
            startIndex += bytesToWrite;
            numBytes -= bytesToWrite;
            assert(numBytes >= 0);
            if (numUncompressedBytes == uncompressedBuffer.length) {
                deflateBlock();
            }
        }
    }

    /**
     * WARNING: flush() affects the output format, because it causes the current contents of uncompressedBuffer
     * to be compressed and written, even if it isn't full.  Unless you know what you're doing, don't call flush().
     * Instead, call close(), which will flush any unwritten data before closing the underlying stream.
     *
     */
    @Override
    public void flush() throws IOException {
        while (numUncompressedBytes > 0) {
            deflateBlock();
        }
        codec.getOutputStream().flush();
    }

    /**
     * close() must be called in order to flush any remaining buffered bytes.
     *
     */
    @Override
    public void close() throws IOException {
        flush();
        if (numberOfThrottleBacks > 0) {
            System.err.println("In BlockCompressedOutputStream, had to throttle back " + numberOfThrottleBacks +
            " times for file " + codec.getOutputFileName());
        }
        codec.close();
    }

    public void write(final int i) throws IOException {
        singleByteArray[0] = (byte)i;
        write(singleByteArray);
    }

    /**
     * Attempt to write the data in uncompressedBuffer to the underlying file in a gzip block.
     * If the entire uncompressedBuffer does not fit in the maximum allowed size, reduce the amount
     * of data to be compressed, and slide the excess down in uncompressedBuffer so it can be picked
     * up in the next deflate event.
     * @return size of gzip block that was written.
     */
    private int deflateBlock() {
        if (numUncompressedBytes == 0) {
            return 0;
        }
        int bytesToCompress = numUncompressedBytes;
        while (true) {
            // Compress the input
            deflater.reset();
            deflater.setInput(uncompressedBuffer, 0, bytesToCompress);
            deflater.finish();
            final int compressedSize = deflater.deflate(compressedBuffer, 0, compressedBuffer.length);

            // If it didn't all fit in compressedBuffer.length, reduce the amount to
            // be compressed and try again.
            if (deflater.getBytesRead() < bytesToCompress) {
                bytesToCompress -= BlockCompressedStreamConstants.UNCOMPRESSED_THROTTLE_AMOUNT;
                ++numberOfThrottleBacks;
                assert(bytesToCompress > 0);
                continue;
            }
            // Data compressed small enough, so write it out.
            crc32.reset();
            crc32.update(uncompressedBuffer, 0, bytesToCompress);
            
            final int totalBlockSize = writeGzipBlock(compressedSize, bytesToCompress, crc32.getValue());
            assert(bytesToCompress <= numUncompressedBytes);

            // Clear out from uncompressedBuffer the data that was written 
            if (bytesToCompress == numUncompressedBytes) {
                numUncompressedBytes = 0;
            } else {
                System.arraycopy(uncompressedBuffer, bytesToCompress, uncompressedBuffer, 0,
                        numUncompressedBytes - bytesToCompress);
                numUncompressedBytes -= bytesToCompress;
            }
            return totalBlockSize;
        }
        // unreachable
    }

    /**
     * Writes the entire gzip block, assuming the compressed data is stored in compressedBuffer
     * @return  size of gzip block that was written.
     */
    private int writeGzipBlock(final int compressedSize, final int uncompressedSize, final long crc) {
        // Init gzip header
        codec.writeByte(BlockCompressedStreamConstants.GZIP_ID1);
        codec.writeByte(BlockCompressedStreamConstants.GZIP_ID2);
        codec.writeByte(BlockCompressedStreamConstants.GZIP_CM_DEFLATE);
        codec.writeByte(BlockCompressedStreamConstants.GZIP_FLG);
        codec.writeInt(0); // Modification time
        codec.writeByte(BlockCompressedStreamConstants.GZIP_XFL);
        codec.writeByte(BlockCompressedStreamConstants.GZIP_OS_UNKNOWN);
        codec.writeShort(BlockCompressedStreamConstants.GZIP_XLEN);
        codec.writeByte(BlockCompressedStreamConstants.BGZF_ID1);
        codec.writeByte(BlockCompressedStreamConstants.BGZF_ID2);
        codec.writeShort(BlockCompressedStreamConstants.BGZF_LEN);
        final int totalBlockSize = compressedSize + BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH +
                BlockCompressedStreamConstants.BLOCK_FOOTER_LENGTH;

        // I don't know why we store block size - 1, but that is what the spec says
        codec.writeShort((short)(totalBlockSize - 1));
        codec.writeBytes(compressedBuffer, 0, compressedSize);
        codec.writeInt((int)crc);
        codec.writeInt(uncompressedSize);
        return totalBlockSize;
    }
}
