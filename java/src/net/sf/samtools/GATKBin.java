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
import java.util.Collections;
import java.util.List;

/**
 * A temporary solution to work around Java access rights issues:
 * override GATKBin and make it public.
 * TODO: Eliminate once we determine the final fate of the BAM index reading code.
 */
public class GATKBin extends Bin {
    public GATKBin(final int referenceSequence, final int binNumber) {
        super(referenceSequence,binNumber);
    }

    public GATKBin(final Bin bin) {
        super(bin.getReferenceSequence(),bin.getBinNumber());
    }

    @Override
    public int getReferenceSequence() {
        return super.getReferenceSequence();
    }

    @Override
    public int getBinNumber() {
        return super.getBinNumber();
    }

    public List<GATKChunk> getGATKChunkList() {
        List<GATKChunk> gatkChunks = new ArrayList<GATKChunk>();
        for(Chunk chunk: getChunkList())
            gatkChunks.add(new GATKChunk(chunk));
        return gatkChunks;
    }

    public void setGATKChunkList(List<GATKChunk> chunks) {
        super.setChunkList(new ArrayList<Chunk>(chunks));
    }
}
