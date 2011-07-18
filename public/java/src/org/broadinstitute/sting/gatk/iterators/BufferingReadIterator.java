/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.LinkedList;
import java.util.NoSuchElementException;
import java.util.Queue;

/**
 * Buffers access to a large stream of reads, replenishing the buffer only when the reads  
 *
 * @author mhanna
 * @version 0.1
 */
public class BufferingReadIterator implements CloseableIterator<SAMRecord> {
    private final CloseableIterator<SAMRecord> wrappedIterator;
    private final Queue<SAMRecord> buffer;
    private final int bufferSize;

    public BufferingReadIterator(final CloseableIterator<SAMRecord> readIterator, final int bufferSize) {
        this.wrappedIterator = readIterator;
        this.buffer = new LinkedList<SAMRecord>();
        this.bufferSize = bufferSize;
    }

    public boolean hasNext() {
        assureBufferFull();
        return !buffer.isEmpty();
    }

    public SAMRecord next() {
        assureBufferFull();
        if(!hasNext()) throw new NoSuchElementException("No next element available");
        return buffer.remove();
    }

    public void close() {
        wrappedIterator.close();
    }

    public void remove() {
        throw new ReviewedStingException("Unable to remove from a BufferingReadIterator");
    }

    /**
     * If the buffer is empty but there are more elements in the iterator,
     */
    private void assureBufferFull() {
        if(!buffer.isEmpty())
            return;
        while(buffer.size() < bufferSize && wrappedIterator.hasNext())
            buffer.add(wrappedIterator.next());
    }
}
