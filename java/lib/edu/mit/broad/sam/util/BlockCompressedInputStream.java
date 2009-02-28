/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2008 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package edu.mit.broad.sam.util;


import java.io.*;
import java.util.zip.GZIPInputStream;

/*
 * Utility class for reading BGZF block compressed files.
 */
public class BlockCompressedInputStream
        extends InputStream
{

    private InputStream mStream = null;
    private RandomAccessFile mFile = null;
    private byte[] mFileBuffer = null;
    private byte[] mCurrentBlock = null;
    private int mCurrentOffset = 0;
    private long mBlockAddress = 0;
    private int mLastBlockLength = 0;


    public BlockCompressedInputStream(final InputStream stream) {
        mStream = toBufferedStream(stream);
        mFile = null;
    }

    public BlockCompressedInputStream(final File file)
        throws IOException {
        mFile = new RandomAccessFile(file, "r");
        mStream = null;
    }

    public int available()
        throws IOException {
        if (mCurrentBlock == null || mCurrentOffset == mCurrentBlock.length) {
            readBlock();
        }
        if (mCurrentBlock == null) {
            return 0;
        }
        return mCurrentBlock.length - mCurrentOffset;
    }

    public void close()
        throws IOException {
        if (mFile != null) {
            mFile.close();
            mFile = null;
        } else if (mStream != null) {
            mStream.close();
            mStream = null;
        }
        // Encourage garbage collection
        mFileBuffer = null;
        mCurrentBlock = null;
    }

    public int read()
        throws IOException {
        return (available() > 0) ? mCurrentBlock[mCurrentOffset++] : -1;
    }

    public int read(final byte[] buffer)
        throws IOException {
        return read(buffer, 0, buffer.length);
    }

    public int read(final byte[] buffer, int offset, int length)
        throws IOException {
        int bytesRead = 0;
        while (length > 0) {
            final int available = available();
            if (available == 0) {
                break;
            }
            final int copyLength = Math.min(length, available);
            System.arraycopy(mCurrentBlock, mCurrentOffset, buffer, offset, copyLength);
            mCurrentOffset += copyLength;
            offset += copyLength;
            length -= copyLength;
            bytesRead += copyLength;
        }
        return bytesRead;
    }

    public void seek(final long pos)
        throws IOException {
        // Note: pos is a special virtual file pointer, not an actual byte offset
        if (mFile == null) {
            throw new IOException("Cannot seek on stream based file");
        }
        // Decode virtual file pointer
        // Upper 48 bits is the byte offset into the compressed stream of a block.
        // Lower 16 bits is the byte offset into the uncompressed stream inside the block.
        final long compressedOffset = pos >> 16;
        final int uncompressedOffset = (int) (pos & 0xFFFF);
        mFile.seek(compressedOffset);
        mBlockAddress = compressedOffset;
        mLastBlockLength = 0;
        readBlock();
        if (uncompressedOffset >= available()) {
            throw new IOException("Invalid file pointer: " + pos);
        }
        mCurrentOffset = uncompressedOffset;
    }

    public long getFilePointer() {
        return ((mBlockAddress << 16) | mCurrentOffset);
    }

    public static boolean isValidFile(final InputStream stream)
        throws IOException {
        if (!stream.markSupported()) {
            throw new RuntimeException("Cannot test non-buffered stream");
        }
        stream.mark(BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH);
        final byte[] buffer = new byte[BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH];
        final int count = readBytes(stream, buffer, 0, BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH);
        stream.reset();
        if (count != BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH) {
            return false;
        }
        return isValidBlockHeader(buffer);
    }

    private static boolean isValidBlockHeader(final byte[] buffer) {
        return (buffer[0] == BlockCompressedStreamConstants.GZIP_ID1 &&
                (buffer[1] & 0xFF) == BlockCompressedStreamConstants.GZIP_ID2 &&
                (buffer[3] & BlockCompressedStreamConstants.GZIP_FLG) != 0 &&
                buffer[10] == BlockCompressedStreamConstants.GZIP_XLEN &&
                buffer[12] == BlockCompressedStreamConstants.BGZF_ID1 &&
                buffer[13] == BlockCompressedStreamConstants.BGZF_ID2);
    }

    private void readBlock()
        throws IOException {

        if (mFileBuffer == null) {
            mFileBuffer = new byte[BlockCompressedStreamConstants.MAX_COMPRESSED_BLOCK_SIZE];
        }
        int count = readBytes(mFileBuffer, 0, BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH);
        if (count == 0) {
            return;
        }
        if (count != BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH) {
            throw new IOException("Premature end of file");
        }
        final int blockLength = unpackInt16(mFileBuffer, BlockCompressedStreamConstants.BLOCK_LENGTH_OFFSET) + 1;
        if (blockLength < BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH || blockLength > mFileBuffer.length) {
            throw new IOException("Unexpected compressed block length: " + blockLength);
        }
        final int remaining = blockLength - BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH;
        count = readBytes(mFileBuffer, BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH, remaining);
        if (count != remaining) {
            throw new IOException("Premature end of file");
        }
        inflateBlock(mFileBuffer, blockLength);
        mCurrentOffset = 0;
        mBlockAddress += mLastBlockLength;
        mLastBlockLength = blockLength; 
    }

    private void inflateBlock(final byte[] compressedBlock, final int compressedLength)
        throws IOException {
        final int uncompressedLength = unpackInt32(compressedBlock, compressedLength-4);
        byte[] buffer = mCurrentBlock;
        mCurrentBlock = null;
        if (buffer == null || buffer.length != uncompressedLength) {
            buffer = new byte[uncompressedLength];
        }
        final GZIPInputStream gzipStream =
            new GZIPInputStream(new ByteArrayInputStream(compressedBlock, 0, compressedLength));
        try {
            final int count = readBytes(gzipStream, buffer, 0, buffer.length);
            if (count != buffer.length) {
                throw new IOException("Block inflate failed");
            }
            // Note: available() does not return zero here.
            // The only safe way to test is to try to read a byte.
            if (gzipStream.read() != -1) {
                throw new IOException("Block inflate failed");
            }
        } finally {
            gzipStream.close();
        }
        mCurrentBlock = buffer;
    }

    private int readBytes(final byte[] buffer, final int offset, final int length)
        throws IOException {
        if (mFile != null) {
            return readBytes(mFile, buffer, offset, length);
        } else if (mStream != null) {
            return readBytes(mStream, buffer, offset, length);
        } else {
            return 0;
        }
    }

    private static int readBytes(final RandomAccessFile file, final byte[] buffer, final int offset, final int length)
        throws IOException {
        int bytesRead = 0;
        while (bytesRead < length) {
            final int count = file.read(buffer, offset + bytesRead, length - bytesRead);
            if (count <= 0) {
                break;
            }
            bytesRead += count;
        }
        return bytesRead;
    }

    private static int readBytes(final InputStream stream, final byte[] buffer, final int offset, final int length)
        throws IOException {
        int bytesRead = 0;
        while (bytesRead < length) {
            final int count = stream.read(buffer, offset + bytesRead, length - bytesRead);
            if (count <= 0) {
                break;
            }
            bytesRead += count;
        }
        return bytesRead;
    }

    private BufferedInputStream toBufferedStream(final InputStream stream) {
        if (stream instanceof BufferedInputStream) {
            return (BufferedInputStream) stream;
        } else {
            return new BufferedInputStream(stream);
        }
    }

    private int unpackInt16(final byte[] buffer, final int offset) {
        return ((buffer[offset] & 0xFF) |
                ((buffer[offset+1] & 0xFF) << 8));
    }

    private int unpackInt32(final byte[] buffer, final int offset) {
        return ((buffer[offset] & 0xFF) |
                ((buffer[offset+1] & 0xFF) << 8) |
                ((buffer[offset+2] & 0xFF) << 16) |
                ((buffer[offset+3] & 0xFF) << 24));
    }
}


