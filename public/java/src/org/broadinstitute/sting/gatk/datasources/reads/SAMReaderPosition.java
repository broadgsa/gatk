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

import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.GATKBAMFileSpan;
import net.sf.samtools.GATKChunk;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.List;

/**
* Created by IntelliJ IDEA.
* User: mhanna
* Date: 10/14/11
* Time: 10:47 PM
* To change this template use File | Settings | File Templates.
*/
class SAMReaderPosition {
    private final SAMReaderID reader;
    private final BlockInputStream inputStream;

    private final List<GATKChunk> positions;
    private PeekableIterator<GATKChunk> positionIterator;

    /**
     * Stores the next block address to read, or -1 if no such block is available.
     */
    private long nextBlockAddress;


    SAMReaderPosition(final SAMReaderID reader, final BlockInputStream inputStream, GATKBAMFileSpan fileSpan) {
        this.reader = reader;
        this.inputStream = inputStream;

        this.positions = fileSpan.getGATKChunks();
        initialize();
    }

    public SAMReaderID getReader() {
        return reader;
    }

    public BlockInputStream getInputStream() {
        return inputStream;
    }

    /**
     * Retrieves the next block address to be read.
     * @return Next block address to be read.
     */
    public long getBlockAddress() {
        return nextBlockAddress;
    }

    /**
     * Retrieves the first offset of interest in the block returned by getBlockAddress().
     * @return First block of interest in this segment.
     */
    public int getFirstOffsetInBlock() {
        return (nextBlockAddress == positionIterator.peek().getBlockStart()) ? positionIterator.peek().getBlockOffsetStart() : 0;
    }

    /**
     * Retrieves the last offset of interest in the block returned by getBlockAddress().
     * @return First block of interest in this segment.
     */
    public int getLastOffsetInBlock() {
        return (nextBlockAddress == positionIterator.peek().getBlockEnd()) ? positionIterator.peek().getBlockOffsetEnd() : 65536;
    }

    public void reset() {
        initialize();
    }

    /**
     * Resets the SAM reader position to its original state.
     */
    private void initialize() {
        this.positionIterator = new PeekableIterator<GATKChunk>(positions.iterator());
        if(positionIterator.hasNext())
            nextBlockAddress = positionIterator.peek().getBlockStart();
        else
            nextBlockAddress = -1;
    }

    /**
     * Advances the current position to the next block to read, given the current position in the file.
     * @param filePosition The current position within the file.
     */
    void advancePosition(final long filePosition) {
        nextBlockAddress = filePosition >> 16;

        // Check the current file position against the iterator; if the iterator is before the current file position,
        // draw the iterator forward.  Remember when performing the check that coordinates are half-open!
        while(positionIterator.hasNext() && isFilePositionPastEndOfChunk(filePosition,positionIterator.peek())) {
            positionIterator.next();

            // If the block iterator has shot past the file pointer, bring the file pointer flush with the start of the current block.
            if(positionIterator.hasNext() && filePosition < positionIterator.peek().getChunkStart()) {
                nextBlockAddress = positionIterator.peek().getBlockStart();
                //System.out.printf("SAMReaderPosition: next block address advanced to %d%n",nextBlockAddress);
                break;
            }
        }

        // If we've shot off the end of the block pointer, notify consumers that iteration is complete.
        if(!positionIterator.hasNext())
            nextBlockAddress = -1;
    }

    private boolean isFilePositionPastEndOfChunk(final long filePosition, final GATKChunk chunk) {
        return filePosition >= chunk.getChunkEnd();
    }
}
