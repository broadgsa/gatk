/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.samtools.GATKBAMFileSpan;
import net.sf.samtools.GATKChunk;
import net.sf.samtools.util.BAMInputStream;
import net.sf.samtools.util.BlockCompressedFilePointerUtil;
import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.RuntimeEOFException;
import net.sf.samtools.util.SeekableStream;
import org.broad.tribble.util.BlockCompressedStreamConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;
import java.util.LinkedList;

/**
 * Presents decompressed blocks to the SAMFileReader.
 */
public class BlockInputStream extends SeekableStream implements BAMInputStream {
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
     * Current position.
     */
    private SAMReaderPosition position;

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
     * Has the buffer been filled since last request?
     */
    private boolean bufferFilled = false;

    /**
     * Create a new block presenting input stream with a dedicated buffer.
     * @param dispatcher the block loading messenger.
     * @param reader the reader for which to load data.
     * @param validate validates the contents read into the buffer against the contents of a Picard BlockCompressedInputStream.
     */
    BlockInputStream(final BGZFBlockLoadingDispatcher dispatcher, final SAMReaderID reader, final boolean validate) {
        this.reader = reader;
        this.length = reader.samFile.length();

        buffer = ByteBuffer.wrap(new byte[64*1024]);
        buffer.order(ByteOrder.LITTLE_ENDIAN);

        // The state of the buffer assumes that the range of data written into the buffer appears in the range
        // [position,limit), while extra capacity exists in the range [limit,capacity)
        buffer.limit(0);

        this.dispatcher = dispatcher;
        // TODO: Kill the region when all we want to do is start at the beginning of the stream and run to the end of the stream.
        this.position = new SAMReaderPosition(reader,this,new GATKBAMFileSpan(new GATKChunk(0,Long.MAX_VALUE)));

        try {
            if(validate) {
                System.out.printf("BlockInputStream %s: BGZF block validation mode activated%n",this);
                validatingInputStream = new BlockCompressedInputStream(reader.samFile);
                // A bug in ValidatingInputStream means that calling getFilePointer() immediately after initialization will result in an NPE.
                // Poke the stream to start reading data.
                validatingInputStream.available();
            }
            else
                validatingInputStream = null;
        }
        catch(IOException ex) {
            throw new ReviewedStingException("Unable to validate against Picard input stream",ex);
        }
    }

    public long length() {
        return length;
    }

    public long getFilePointer() {
        long filePointer;
        synchronized(lock) {
            if(buffer.remaining() > 0) {
                // If there's data in the buffer, figure out from whence it came.
                final long blockAddress = blockPositions.size() > 0 ? blockPositions.get(0) : 0;
                final int blockOffset = buffer.position();
                filePointer = blockAddress << 16 | blockOffset;
            }
            else {
                // Otherwise, find the next position to load.
                filePointer = position.getBlockAddress() << 16;
            }
        }

        if(validatingInputStream != null && filePointer != validatingInputStream.getFilePointer())
            throw new ReviewedStingException(String.format("Position of input stream is invalid; expected (block address, block offset) = (%d,%d), got (%d,%d)",
                    BlockCompressedFilePointerUtil.getBlockAddress(filePointer),BlockCompressedFilePointerUtil.getBlockOffset(filePointer),
                    BlockCompressedFilePointerUtil.getBlockAddress(validatingInputStream.getFilePointer()),BlockCompressedFilePointerUtil.getBlockOffset(validatingInputStream.getFilePointer())));

        return filePointer;
    }

    public void seek(long target) {
        // TODO: Validate the seek point.
        //System.out.printf("Thread %s, BlockInputStream %s: seeking to block %d, offset %d%n",Thread.currentThread().getId(),this,BlockCompressedFilePointerUtil.getBlockAddress(target),BlockCompressedFilePointerUtil.getBlockOffset(target));
        synchronized(lock) {
            clearBuffers();
            position.advancePosition(BlockCompressedFilePointerUtil.getBlockAddress(target));
            waitForBufferFill();
            buffer.position(BlockCompressedFilePointerUtil.getBlockOffset(target));

            if(validatingInputStream != null) {
                try {
                    validatingInputStream.seek(target);
                }
                catch(IOException ex) {
                    throw new ReviewedStingException("Unable to validate against Picard input stream",ex);
                }
            }
        }
    }

    private void clearBuffers() {
        this.position.reset();

        // Buffer semantics say that outside of a lock, buffer should always be prepared for reading.
        // Indicate no data to be read.
        buffer.clear();
        buffer.limit(0);

        blockOffsets.clear();
        blockPositions.clear();
    }

    public boolean eof() {
        synchronized(lock) {
            // TODO: Handle multiple empty BGZF blocks at end of the file.
            return position != null && position.getBlockAddress() >= length;
        }
    }

    public void setCheckCrcs(final boolean check) {
        // TODO: Implement
    }

    /**
     * Submits a new access plan for the given dataset.
     * @param position The next seek point for BAM data in this reader.
     */
    public void submitAccessPlan(final SAMReaderPosition position) {
        //System.out.printf("Thread %s: submitting access plan for block at position: %d%n",Thread.currentThread().getId(),position.getBlockAddress());
        synchronized(lock) {
            // Assume that the access plan is going to tell us to start where we are and move forward.
            // If this isn't the case, we'll soon receive a seek request and the buffer will be forced to reset.
            if(this.position != null && position.getBlockAddress() < this.position.getBlockAddress())
                position.advancePosition(this.position.getBlockAddress());
        }
        this.position = position;
    }

    private void compactBuffer() {
        // Compact buffer to maximize storage space.
        int bytesToRemove = 0;

        // Look ahead to see if we can compact away the first block in the series.
        while(blockOffsets.size() > 1 && buffer.position() < blockOffsets.get(1)) {
            bytesToRemove += blockOffsets.remove();
            blockPositions.remove();
        }

        // If we end up with an empty block at the end of the series, compact this as well.
        if(buffer.remaining() == 0 && !blockOffsets.isEmpty() && buffer.position() >= blockOffsets.peek()) {
            bytesToRemove += buffer.position();
            blockOffsets.remove();
            blockPositions.remove();
        }

        int finalBufferStart = buffer.position() - bytesToRemove;
        int finalBufferSize = buffer.remaining();

        buffer.position(bytesToRemove);
        buffer.compact();

        buffer.position(finalBufferStart);
        buffer.limit(finalBufferStart+finalBufferSize);
    }

    /**
     * Push contents of incomingBuffer into the end of this buffer.
     * MUST be called from a thread that is NOT the reader thread.
     * @param incomingBuffer The data being pushed into this input stream.
     * @param position target position for the data.
     */
    public void copyIntoBuffer(final ByteBuffer incomingBuffer, final SAMReaderPosition position, final long filePosition) {
        synchronized(lock) {
            try {
                compactBuffer();
                // Open up the buffer for more reading.
                buffer.limit(buffer.capacity());

                // Advance the position to take the most recent read into account.
                long lastReadPosition = position.getBlockAddress();

                byte[] validBytes = null;
                if(validatingInputStream != null) {
                    validBytes = new byte[incomingBuffer.remaining()];

                    byte[] currentBytes = new byte[incomingBuffer.remaining()];
                    int pos = incomingBuffer.position();
                    int lim = incomingBuffer.limit();
                    incomingBuffer.get(currentBytes);

                    incomingBuffer.limit(lim);
                    incomingBuffer.position(pos);

                    long currentFilePointer = validatingInputStream.getFilePointer();
                    validatingInputStream.seek(lastReadPosition << 16);
                    validatingInputStream.read(validBytes);
                    validatingInputStream.seek(currentFilePointer);

                    if(!Arrays.equals(validBytes,currentBytes))
                        throw new ReviewedStingException(String.format("Bytes being inserted into BlockInputStream %s are incorrect",this));
                }

                this.position = position;
                position.advancePosition(filePosition);

                if(buffer.remaining() < incomingBuffer.remaining()) {
                    //System.out.printf("Thread %s: waiting for available space in buffer; buffer remaining = %d, incoming buffer remaining = %d%n",Thread.currentThread().getId(),buffer.remaining(),incomingBuffer.remaining());
                    lock.wait();
                    //System.out.printf("Thread %s: waited for available space in buffer; buffer remaining = %d, incoming buffer remaining = %d%n", Thread.currentThread().getId(), buffer.remaining(), incomingBuffer.remaining());
                }

                // Queue list of block offsets / block positions.
                blockOffsets.add(buffer.position());
                blockPositions.add(lastReadPosition);

                buffer.put(incomingBuffer);

                // Set up the buffer for reading.
                buffer.flip();
                bufferFilled = true;

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
                ReviewedStingException toThrow = new ReviewedStingException(String.format("Thread %s, BlockInputStream %s: Unable to retrieve BAM data from disk",Thread.currentThread().getId(),this),error);
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

        if(validatingInputStream != null) {
            byte[] validBytes = new byte[length];
            try {
                validatingInputStream.read(validBytes,offset,length);
                for(int i = offset; i < offset+length; i++) {
                    if(bytes[i] != validBytes[i]) {
                        System.out.printf("Thread %s: preparing to throw an exception because contents don't match%n",Thread.currentThread().getId());
                        throw new ReviewedStingException(String.format("Thread %s: blockInputStream %s attempting to return wrong set of bytes; mismatch at offset %d",Thread.currentThread().getId(),this,i));
                    }
                }
            }
            catch(IOException ex) {
                throw new ReviewedStingException("Unable to validate against Picard input stream",ex);
            }
        }

        return length - remaining;
    }

    public void close() {
        if(validatingInputStream != null) {
            try {
                validatingInputStream.close();
            }
            catch(IOException ex) {
                throw new ReviewedStingException("Unable to validate against Picard input stream",ex);
            }
        }
    }

    public String getSource() {
        return reader.getSamFilePath();
    }

    private void waitForBufferFill() {
        synchronized(lock) {
            bufferFilled = false;
            if(buffer.remaining() == 0 && !eof()) {
                //System.out.printf("Thread %s is waiting for a buffer fill from position %d to buffer %s%n",Thread.currentThread().getId(),position.getBlockAddress(),this);
                dispatcher.queueBlockLoad(position);
                try {
                    lock.wait();
                }
                catch(InterruptedException ex) {
                    // TODO: handle me.
                    throw new ReviewedStingException("Interrupt occurred waiting for buffer to fill",ex);
                }

                if(bufferFilled && buffer.remaining() == 0)
                    throw new RuntimeEOFException("No more data left in InputStream");
            }
        }
    }
}
