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

import htsjdk.samtools.GATKBAMFileSpan;
import htsjdk.samtools.GATKChunk;
import htsjdk.samtools.util.BlockCompressedInputStream;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;

import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * Presents decompressed blocks to the SAMFileReader.
 */
public class BlockInputStream extends InputStream {
    /**
     * Mechanism for triggering block loads.
     */
    private final BGZFBlockLoadingDispatcher dispatcher;

    /**
     * The reader whose data is supplied by this input stream.
     */
    private final SAMReaderID reader;

    /**
     * Length of the input stream.
     */
    private final long length;

    /**
     * The latest error reported by an asynchronous block load.
     */
    private Throwable error;

    /**
     * Current accessPlan.
     */
    private BAMAccessPlan accessPlan;

    /**
     * A stream of compressed data blocks.
     */
    private final ByteBuffer buffer;

    /**
     * Offsets of the given blocks in the buffer.
     */
    private LinkedList<Integer> blockOffsets = new LinkedList<Integer>();

    /**
     * Source positions of the given blocks in the buffer.
     */
    private LinkedList<Long> blockPositions = new LinkedList<Long>();

    /**
     * Provides a lock to wait for more data to arrive.
     */
    private final Object lock = new Object();

    /**
     * An input stream to use when comparing data back to what it should look like.
     */
    private final BlockCompressedInputStream validatingInputStream;

    /**
     * Create a new block presenting input stream with a dedicated buffer.
     * @param dispatcher the block loading messenger.
     * @param reader the reader for which to load data.
     * @param validate validates the contents read into the buffer against the contents of a Picard BlockCompressedInputStream.
     */
    BlockInputStream(final BGZFBlockLoadingDispatcher dispatcher, final SAMReaderID reader, final boolean validate) {
        this.reader = reader;
        this.length = reader.getSamFile().length();

        buffer = ByteBuffer.wrap(new byte[64*1024]);
        buffer.order(ByteOrder.LITTLE_ENDIAN);

        // The state of the buffer assumes that the range of data written into the buffer appears in the range
        // [position,limit), while extra capacity exists in the range [limit,capacity)
        buffer.limit(0);

        this.dispatcher = dispatcher;
        // TODO: Kill the region when all we want to do is start at the beginning of the stream and run to the end of the stream.
        this.accessPlan = new BAMAccessPlan(reader,this,new GATKBAMFileSpan(new GATKChunk(0,Long.MAX_VALUE)));

        // The block offsets / block positions guarantee that the ending offset/position in the data structure maps to
        // the point in the file just following the last read.  These two arrays should never be empty; initializing
        // to 0 to match the position above.
        this.blockOffsets.add(0);
        this.blockPositions.add(0L);

        try {
            if(validate) {
                System.out.printf("BlockInputStream %s: BGZF block validation mode activated%n",this);
                validatingInputStream = new BlockCompressedInputStream(reader.getSamFile());
                // A bug in ValidatingInputStream means that calling getFilePointer() immediately after initialization will result in an NPE.
                // Poke the stream to start reading data.
                validatingInputStream.available();
            }
            else
                validatingInputStream = null;
        }
        catch(IOException ex) {
            throw new ReviewedGATKException("Unable to validate against Picard input stream",ex);
        }
    }

    public long length() {
        return length;
    }

    public long getFilePointer() {
        long filePointer;
        synchronized(lock) {
            // Find the current block within the input stream.
            int blockIndex;
            for(blockIndex = 0; blockIndex+1 < blockOffsets.size() && buffer.position() > blockOffsets.get(blockIndex+1); blockIndex++)
                ;
            filePointer = blockPositions.get(blockIndex) + (buffer.position()-blockOffsets.get(blockIndex));
        }

//        if(validatingInputStream != null && filePointer != validatingInputStream.getFilePointer())
//            throw new ReviewedGATKException(String.format("Position of input stream is invalid; expected (block address, block offset) = (%d,%d), got (%d,%d)",
//                    BlockCompressedFilePointerUtil.getBlockAddress(validatingInputStream.getFilePointer()),BlockCompressedFilePointerUtil.getBlockOffset(validatingInputStream.getFilePointer()),
//                    BlockCompressedFilePointerUtil.getBlockAddress(filePointer),BlockCompressedFilePointerUtil.getBlockOffset(filePointer)));

        return filePointer;
    }

    private void clearBuffers() {
        this.accessPlan.reset();

        // Buffer semantics say that outside of a lock, buffer should always be prepared for reading.
        // Indicate no data to be read.
        buffer.clear();
        buffer.limit(0);

        // Clear everything except the last block offset / position
        blockOffsets.clear();
        blockOffsets.add(0);
        while(blockPositions.size() > 1)
            blockPositions.removeFirst();
    }

    public boolean eof() {
        synchronized(lock) {
            // TODO: Handle multiple empty BGZF blocks at end of the file.
            return accessPlan != null && (accessPlan.getBlockAddress() < 0 || accessPlan.getBlockAddress() >= length);
        }
    }

    /**
     * Submits a new access plan for the given dataset and seeks to the given point.
     * @param accessPlan The next seek point for BAM data in this reader.
     */
    public void submitAccessPlan(final BAMAccessPlan accessPlan) {
        //System.out.printf("Thread %s: submitting access plan for block at position: %d%n",Thread.currentThread().getId(),position.getBlockAddress());
        this.accessPlan = accessPlan;
        accessPlan.reset();

        clearBuffers();

        // Pull the iterator past any oddball chunks at the beginning of the shard (chunkEnd < chunkStart, empty chunks, etc).
        // TODO: Don't pass these empty chunks in.
        accessPlan.advancePosition(makeFilePointer(accessPlan.getBlockAddress(),0));

        if(accessPlan.getBlockAddress() >= 0) {
            waitForBufferFill();
        }

        if(validatingInputStream != null) {
            try {
                validatingInputStream.seek(makeFilePointer(accessPlan.getBlockAddress(),0));
            }
            catch(IOException ex) {
                throw new ReviewedGATKException("Unable to validate against Picard input stream",ex);
            }
        }

    }


    private void compactBuffer() {
        // Compact buffer to maximize storage space.
        int bytesToRemove = 0;

        // Look ahead to see if we can compact away the first blocks in the series.
        while(blockOffsets.size() > 1 && buffer.position() >= blockOffsets.get(1)) {
            blockOffsets.remove();
            blockPositions.remove();
            bytesToRemove = blockOffsets.peek();
        }

        // If we end up with an empty block at the end of the series, compact this as well.
        if(buffer.remaining() == 0 && blockOffsets.size() > 1 && buffer.position() >= blockOffsets.peek()) {
            bytesToRemove += buffer.position();
            blockOffsets.remove();
            blockPositions.remove();
        }

        int finalBufferStart = buffer.position() - bytesToRemove;
        int finalBufferSize = buffer.remaining();

        // Position the buffer to remove the unneeded data, and compact it away.
        buffer.position(bytesToRemove);
        buffer.compact();

        // Reset the limits for reading.
        buffer.position(finalBufferStart);
        buffer.limit(finalBufferStart+finalBufferSize);

        // Shift everything in the offset buffer down to accommodate the bytes removed from the buffer.
        for(int i = 0; i < blockOffsets.size(); i++)
            blockOffsets.set(i,blockOffsets.get(i)-bytesToRemove);
    }

    /**
     * Push contents of incomingBuffer into the end of this buffer.
     * MUST be called from a thread that is NOT the reader thread.
     * @param incomingBuffer The data being pushed into this input stream.
     * @param accessPlan target access plan for the data.
     * @param filePosition the current position of the file pointer
     */
    public void copyIntoBuffer(final ByteBuffer incomingBuffer, final BAMAccessPlan accessPlan, final long filePosition) {
        synchronized(lock) {
            try {
                if(validatingInputStream != null) {
                    byte[] validBytes = new byte[incomingBuffer.remaining()];

                    byte[] currentBytes = new byte[incomingBuffer.remaining()];
                    int pos = incomingBuffer.position();
                    int lim = incomingBuffer.limit();
                    incomingBuffer.get(currentBytes);

                    incomingBuffer.limit(lim);
                    incomingBuffer.position(pos);

                    long currentFilePointer = validatingInputStream.getFilePointer();
                    validatingInputStream.seek(makeFilePointer(accessPlan.getBlockAddress(), 0));
                    validatingInputStream.read(validBytes);
                    validatingInputStream.seek(currentFilePointer);

                    if(!Arrays.equals(validBytes,currentBytes))
                        throw new ReviewedGATKException(String.format("Bytes being inserted into BlockInputStream %s are incorrect",this));
                }

                compactBuffer();
                // Open up the buffer for more reading.
                buffer.limit(buffer.capacity());

                // Get the spans overlapping this particular block...
                List<GATKChunk> spansOverlapping = accessPlan.getSpansOverlappingBlock(accessPlan.getBlockAddress(),filePosition);

                // ...and advance the block
                this.accessPlan = accessPlan;
                accessPlan.advancePosition(makeFilePointer(filePosition, 0));

                if(buffer.remaining() < incomingBuffer.remaining())
                    lock.wait();

                final int bytesInIncomingBuffer = incomingBuffer.limit();

                for(GATKChunk spanOverlapping: spansOverlapping) {
                    // Clear out the endcap tracking state and add in the starting position for this transfer.
                    blockOffsets.removeLast();
                    blockOffsets.add(buffer.position());
                    blockPositions.removeLast();
                    blockPositions.add(spanOverlapping.getChunkStart());

                    // Stream the buffer into the data stream.
                    incomingBuffer.limit((spanOverlapping.getBlockEnd() > spanOverlapping.getBlockStart()) ? bytesInIncomingBuffer : spanOverlapping.getBlockOffsetEnd());
                    incomingBuffer.position(spanOverlapping.getBlockOffsetStart());
                    buffer.put(incomingBuffer);

                    // Add the endcap for this transfer.
                    blockOffsets.add(buffer.position());
                    blockPositions.add(spanOverlapping.getChunkEnd());
                }

                // Set up the buffer for reading.
                buffer.flip();

                lock.notify();
            }
            catch(Exception ex) {
                reportException(ex);
                lock.notify();
            }
        }
    }

    void reportException(Throwable t) {
        synchronized(lock) {
            this.error = t;
            lock.notify();
        }
    }

    private void checkForErrors() {
        synchronized(lock) {
            if(error != null) {
                ReviewedGATKException toThrow = new ReviewedGATKException(String.format("Thread %s, BlockInputStream %s: Unable to retrieve BAM data from disk",Thread.currentThread().getId(),this),error);
                toThrow.setStackTrace(error.getStackTrace());
                throw toThrow;
            }
        }
    }

    /**
     * Reads the next byte of data from the input stream.
     * @return Next byte of data, from 0->255, as an int.
     */
    @Override
    public int read() {
        byte[] singleByte = new byte[1];
        read(singleByte);
        return singleByte[0];
    }

    /**
     * Fills the given byte array to the extent possible.
     * @param bytes byte array to be filled.
     * @return The number of bytes actually read.
     */
    @Override
    public int read(byte[] bytes) {
        return read(bytes,0,bytes.length);
    }

    @Override
    public int read(byte[] bytes, final int offset, final int length) {
        int remaining = length;
        synchronized(lock) {
            while(remaining > 0) {
                // Check for error conditions during last read.
                checkForErrors();

                // If completely out of space, queue up another buffer fill.
                waitForBufferFill();

                // Couldn't manage to load any data at all; abort and return what's available.
                if(buffer.remaining() == 0)
                    break;

                int numBytesToCopy = Math.min(buffer.remaining(),remaining);
                buffer.get(bytes,length-remaining+offset,numBytesToCopy);
                remaining -= numBytesToCopy;

                //if(remaining > 0)
                //    System.out.printf("Thread %s: read the first %d bytes of a %d byte request%n",Thread.currentThread().getId(),length-remaining,length);
                // TODO: Assert that we don't copy across a block boundary
            }

            // Notify any waiting threads that some of the contents of the buffer were removed.
            if(length-remaining > 0)
                lock.notify();
        }

//        if(validatingInputStream != null) {
//            byte[] validBytes = new byte[length];
//            try {
//                validatingInputStream.read(validBytes,offset,length);
//                for(int i = offset; i < offset+length; i++) {
//                    if(bytes[i] != validBytes[i])
//                        throw new ReviewedGATKException(String.format("Thread %s: blockInputStream %s attempting to return wrong set of bytes; mismatch at offset %d",Thread.currentThread().getId(),this,i));
//                }
//            }
//            catch(IOException ex) {
//                throw new ReviewedGATKException("Unable to validate against Picard input stream",ex);
//            }
//        }

        // If any data was copied into the buffer, return the amount of data copied.
        if(remaining < length)
            return length - remaining;

        // Otherwise, return -1.
        return -1;
    }

    public void close() {
        if(validatingInputStream != null) {
            try {
                validatingInputStream.close();
            }
            catch(IOException ex) {
                throw new ReviewedGATKException("Unable to validate against Picard input stream",ex);
            }
        }
    }

    public String getSource() {
        return reader.getSamFilePath();
    }

    private void waitForBufferFill() {
        synchronized(lock) {
            if(buffer.remaining() == 0 && !eof()) {
                //System.out.printf("Thread %s is waiting for a buffer fill from position %d to buffer %s%n",Thread.currentThread().getId(),position.getBlockAddress(),this);
                dispatcher.queueBlockLoad(accessPlan);
                try {
                    lock.wait();
                }
                catch(InterruptedException ex) {
                    throw new ReviewedGATKException("Interrupt occurred waiting for buffer to fill",ex);
                }
            }
        }
    }

    /**
     * Create an encoded BAM file pointer given the address of a BGZF block and an offset.
     * @param blockAddress Physical address on disk of a BGZF block.
     * @param blockOffset Offset into the uncompressed data stored in the BGZF block.
     * @return 64-bit pointer encoded according to the BAM spec.
     */
    public static long makeFilePointer(final long blockAddress, final int blockOffset) {
        return blockAddress << 16 | blockOffset;
    }
}
