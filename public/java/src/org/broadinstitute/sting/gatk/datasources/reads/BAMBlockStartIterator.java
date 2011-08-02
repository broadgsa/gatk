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

import org.broadinstitute.sting.utils.exceptions.StingException;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Feb 7, 2011
 * Time: 2:46:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class BAMBlockStartIterator implements Iterator<Long> {
    /**
     * How large is a BGZF header?
     */
    private static int BGZF_HEADER_SIZE = 18;

    /**
     * Where within the header does the BLOCKSIZE actually live?
     */
    private static int BLOCK_SIZE_HEADER_POSITION = BGZF_HEADER_SIZE - 2;

    private FileChannel bamInputChannel;
    private ByteBuffer headerByteBuffer;

    private long nextLocation = 0;

    public BAMBlockStartIterator(File bamFile) {
        try {
            FileInputStream bamInputStream = new FileInputStream(bamFile);
            bamInputChannel = bamInputStream.getChannel();

            headerByteBuffer = ByteBuffer.allocate(BGZF_HEADER_SIZE);
            headerByteBuffer.order(ByteOrder.LITTLE_ENDIAN);

        }
        catch(IOException ex) {
            throw new StingException("Could not open file",ex);
        }
    }

    public boolean hasNext() {
        return nextLocation != -1;
    }

    public Long next() {
        long currentLocation = nextLocation;
        advance();
        return currentLocation;
    }

    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a BAMBlockStartIterator");
    }

    private void advance() {
        int readStatus;

        headerByteBuffer.clear();
        try {
            readStatus = bamInputChannel.read(headerByteBuffer);
        }
        catch(IOException ex) {
            throw new StingException("Could not read header data",ex);
        }

        if(readStatus == -1) {
            nextLocation = -1;
            try {
                bamInputChannel.close();
            }
            catch(IOException ex) {
                throw new StingException("Could not close input file",ex);
            }
            return;
        }

        headerByteBuffer.position(BLOCK_SIZE_HEADER_POSITION);
        int blockSize = headerByteBuffer.getShort();

        try {
            bamInputChannel.position(bamInputChannel.position()+blockSize-BGZF_HEADER_SIZE+1);
            nextLocation = bamInputChannel.position();
        }
        catch(IOException ex) {
            throw new StingException("Could not reposition input stream",ex);
        }
    }

    public static void main(String argv[]) throws IOException {
        BAMBlockStartIterator blockStartIterator = new BAMBlockStartIterator(new File("/Users/mhanna/testdata/reads/MV1994.bam"));
        int i = 0;
        while(blockStartIterator.hasNext())
            System.out.printf("%d -> %d%n",i++,blockStartIterator.next());
    }
}
