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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A temporary solution to work around Java access rights issues:
 * override BAMFileSpan and make it public.
 * TODO: Eliminate once we determine the final fate of the BAM index reading code.
 */
public class GATKBAMFileSpan extends BAMFileSpan {
    /**
     * Create a new empty list of chunks.
     */
    public GATKBAMFileSpan() {
        super();
    }

    /**
     * Convenience constructor to construct a BAM file span from
     * a single chunk.
     * @param chunk Chunk to use as the sole region in this span.
     */
    public GATKBAMFileSpan(final Chunk chunk) {
        super(chunk);
    }

    /**
     * Create a new chunk list from the given list of chunks.
     * @param chunks Constituent chunks.
     */
    public GATKBAMFileSpan(final GATKChunk[] chunks) {
        super(Arrays.<Chunk>asList(chunks));
    }

    /**
     * Gets the constituent chunks stored in this span.
     * @return An unmodifiable list of chunks.
     */
    public List<GATKChunk> getGATKChunks() {
        List<GATKChunk> gatkChunks = new ArrayList<GATKChunk>();
        for(Chunk chunk: getChunks())
            gatkChunks.add(new GATKChunk(chunk));
        return gatkChunks;
    }

    public String toString() {
        StringBuilder builder = new StringBuilder();
        for(GATKChunk chunk: getGATKChunks())
            builder.append(String.format("%s;",chunk));
        return builder.toString();
    }
}
