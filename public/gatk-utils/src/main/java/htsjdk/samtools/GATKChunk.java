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

package htsjdk.samtools;

/**
 * A temporary solution to work around Java access rights issues:
 * override chunk and make it public.
 * TODO: Eliminate once we determine the final fate of the BAM index reading code.
 */
public class GATKChunk extends Chunk {
    /**
     * The average ratio of compressed block size / uncompressed block size, computed empirically
     * using the output of org.broadinstitute.gatk.engine.datasources.reads.utilities.PrintBGZFBounds.
     */
    private static final double AVERAGE_BAM_COMPRESSION_RATIO = 0.39;

    public GATKChunk(final long start, final long stop) {
        super(start,stop);
    }

    public GATKChunk(final long blockStart, final int blockOffsetStart, final long blockEnd, final int blockOffsetEnd) {
        super(blockStart << 16 | blockOffsetStart,blockEnd << 16 | blockOffsetEnd);
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

    public long getBlockStart() {
        return getChunkStart() >>> 16;
    }

    public int getBlockOffsetStart() {
        return (int)(getChunkStart() & 0xFFFF);
    }

    public long getBlockEnd() {
        return getChunkEnd() >>> 16;
    }

    public int getBlockOffsetEnd() {
        return ((int)getChunkEnd() & 0xFFFF);
    }

    /**
     * Computes an approximation of the uncompressed size of the
     * chunk, in bytes.  Can be used to determine relative weights
     * of chunk size.
     * @return An approximation of the chunk size in bytes.
     */
    public long size() {
        final long chunkSpan = Math.round(((getChunkEnd()>>16)-(getChunkStart()>>16))/AVERAGE_BAM_COMPRESSION_RATIO);
        final int offsetSpan = (int)((getChunkEnd()&0xFFFF)-(getChunkStart()&0xFFFF));
        return chunkSpan + offsetSpan;
    }

    /**
     * Merges two chunks together. The caller is responsible for testing whether the
     * chunks overlap/are adjacent before calling this method!
     *
     * @param other the chunk to merge with this chunk
     * @return a new chunk representing the union of the two chunks (provided the chunks were
     *         overlapping/adjacent)
     */
    public GATKChunk merge ( GATKChunk other ) {
        return new GATKChunk(Math.min(getChunkStart(), other.getChunkStart()), Math.max(getChunkEnd(), other.getChunkEnd()));
    }
}
