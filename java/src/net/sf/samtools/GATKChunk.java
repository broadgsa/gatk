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

package net.sf.samtools;

/**
 * A temporary solution to work around Java access rights issues:
 * override chunk and make it public.
 * TODO: Eliminate once we determine the final fate of the BAM index reading code.
 */
public class GATKChunk extends Chunk {
    public GATKChunk(final long start, final long stop) {
        super(start,stop);
    }

    public GATKChunk(final Chunk chunk) {
        super(chunk.getChunkStart(),chunk.getChunkEnd());
    }

    @Override
    public GATKChunk clone() {
        return new GATKChunk(getChunkStart(),getChunkEnd());
    }

    @Override
    public long getChunkStart() {
        return super.getChunkStart();
    }

    @Override
    public void setChunkStart(final long value) {
        super.setChunkStart(value);
    }

    @Override
    public long getChunkEnd() {
        return super.getChunkEnd();
    }

    @Override
    public void setChunkEnd(final long value) {
        super.setChunkEnd(value);
    }


}
