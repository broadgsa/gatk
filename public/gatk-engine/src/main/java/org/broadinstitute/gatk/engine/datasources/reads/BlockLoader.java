/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.datasources.reads;

import htsjdk.samtools.util.BlockCompressedStreamConstants;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * An engine for loading blocks.
 */
class BlockLoader implements Runnable {
    /**
     * Coordinates the input queue.
     */
    private BGZFBlockLoadingDispatcher dispatcher;

    /**
     * A cache from which to retrieve open file handles.
     */
    private final FileHandleCache fileHandleCache;

    /**
     * Whether asynchronous decompression should happen.
     */
    private final boolean decompress;

    /**
     * An direct input buffer for incoming data from disk.
     */
    private final ByteBuffer inputBuffer;

    public BlockLoader(final BGZFBlockLoadingDispatcher dispatcher, final FileHandleCache fileHandleCache, final boolean decompress) {
        this.dispatcher = dispatcher;
        this.fileHandleCache = fileHandleCache;
        this.decompress = decompress;

        this.inputBuffer = ByteBuffer.allocateDirect(64*1024 + BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length);
        inputBuffer.order(ByteOrder.LITTLE_ENDIAN);
    }

    public void run() {
        for(;;) {
            BAMAccessPlan accessPlan = null;
            try {
                accessPlan = dispatcher.claimNextWorkRequest();
                FileInputStream inputStream = fileHandleCache.claimFileInputStream(accessPlan.getReader());

                //long blockAddress = readerPosition.getBlockAddress();
                //System.out.printf("Thread %s: BlockLoader: copying bytes from %s at position %d into %s%n",Thread.currentThread().getId(),inputStream,blockAddress,readerPosition.getInputStream());

                ByteBuffer compressedBlock = readBGZFBlock(inputStream,accessPlan.getBlockAddress());
                long nextBlockAddress = position(inputStream);
                fileHandleCache.releaseFileInputStream(accessPlan.getReader(),inputStream);

                ByteBuffer block = decompress ? decompressBGZFBlock(compressedBlock) : compressedBlock;
                int bytesCopied = block.remaining();

                BlockInputStream bamInputStream = accessPlan.getInputStream();
                bamInputStream.copyIntoBuffer(block,accessPlan,nextBlockAddress);

                //System.out.printf("Thread %s: BlockLoader: copied %d bytes from %s at position %d into %s%n",Thread.currentThread().getId(),bytesCopied,inputStream,blockAddress,readerPosition.getInputStream());
            }
            catch(Throwable error) {
                if(accessPlan != null && accessPlan.getInputStream() != null)
                    accessPlan.getInputStream().reportException(error);
            }
        }

    }

    private ByteBuffer readBGZFBlock(final FileInputStream inputStream, final long blockAddress) throws IOException {
        FileChannel channel = inputStream.getChannel();

        // Read the block header
        channel.position(blockAddress);

        int uncompressedDataSize = 0;
        int bufferSize = 0;

        do {
            inputBuffer.clear();
            inputBuffer.limit(BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH);
            channel.read(inputBuffer);

            // Read out the size of the full BGZF block into a two bit short container, then 'or' that
            // value into an int buffer to transfer the bitwise contents into an int.
            inputBuffer.flip();
            if(inputBuffer.remaining() != BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH)
                throw new ReviewedGATKException("BUG: unable to read a the complete block header in one pass.");

            // Verify that the file was read at a valid point.
            if(unpackUByte8(inputBuffer,0) != BlockCompressedStreamConstants.GZIP_ID1 ||
                    unpackUByte8(inputBuffer,1) != BlockCompressedStreamConstants.GZIP_ID2 ||
                    unpackUByte8(inputBuffer,3) != BlockCompressedStreamConstants.GZIP_FLG ||
                    unpackUInt16(inputBuffer,10) != BlockCompressedStreamConstants.GZIP_XLEN ||
                    unpackUByte8(inputBuffer,12) != BlockCompressedStreamConstants.BGZF_ID1 ||
                    unpackUByte8(inputBuffer,13) != BlockCompressedStreamConstants.BGZF_ID2) {
                throw new ReviewedGATKException("BUG: Started reading compressed block at incorrect position");
            }

            inputBuffer.position(BlockCompressedStreamConstants.BLOCK_LENGTH_OFFSET);
            bufferSize = unpackUInt16(inputBuffer,BlockCompressedStreamConstants.BLOCK_LENGTH_OFFSET)+1;

            // Adjust buffer limits and finish reading the block.  Also read the next header, just in case there's a 0-byte block.
            inputBuffer.limit(bufferSize);
            inputBuffer.position(BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH);
            channel.read(inputBuffer);

            // Check the uncompressed length.  If 0 and not at EOF, we'll want to check the next block.
            uncompressedDataSize = inputBuffer.getInt(inputBuffer.limit()-4);
            //System.out.printf("Uncompressed block size of the current block (at position %d) is %d%n",channel.position()-inputBuffer.limit(),uncompressedDataSize);
        }
        while(uncompressedDataSize == 0 && channel.position() < channel.size());

        // Prepare the buffer for reading.
        inputBuffer.flip();

        return inputBuffer;
    }

    private ByteBuffer decompressBGZFBlock(final ByteBuffer bgzfBlock) throws DataFormatException {
        final int compressedBufferSize = bgzfBlock.remaining();

        // Determine the uncompressed buffer size (
        bgzfBlock.position(bgzfBlock.limit()-4);
        int uncompressedBufferSize = bgzfBlock.getInt();
        byte[] uncompressedContent = new byte[uncompressedBufferSize];

        // Bound the CDATA section of the buffer.
        bgzfBlock.limit(compressedBufferSize-BlockCompressedStreamConstants.BLOCK_FOOTER_LENGTH);
        bgzfBlock.position(BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH);
        byte[] compressedContent = new byte[bgzfBlock.remaining()];
        ByteBuffer.wrap(compressedContent).put(bgzfBlock);

        // Decompress the buffer.
        final Inflater inflater = new Inflater(true);
        inflater.setInput(compressedContent);
        int bytesUncompressed = inflater.inflate(uncompressedContent);
        if(bytesUncompressed != uncompressedBufferSize)
            throw new ReviewedGATKException("Error decompressing block");

        return ByteBuffer.wrap(uncompressedContent);
    }

    private long position(final FileInputStream inputStream) throws IOException {
        return inputStream.getChannel().position();
    }

    private int unpackUByte8(final ByteBuffer buffer,final int position) {
        return buffer.get(position) & 0xFF;
    }

    private int unpackUInt16(final ByteBuffer buffer,final int position) {
        // Read out the size of the full BGZF block into a two bit short container, then 'or' that
        // value into an int buffer to transfer the bitwise contents into an int.
        return buffer.getShort(position) & 0xFFFF;
    }
}
