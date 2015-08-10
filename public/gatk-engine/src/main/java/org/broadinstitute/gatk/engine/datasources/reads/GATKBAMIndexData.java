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

import htsjdk.samtools.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Stores and processes a single reference worth of GATK data.
 */
public class GATKBAMIndexData {
    private final GATKBAMIndex index;
    private final int referenceSequence;
    private final List<GATKBin> bins;
    private final LinearIndex linearIndex;

    public GATKBAMIndexData(final GATKBAMIndex index, final int referenceSequence, final List<GATKBin> bins, final LinearIndex linearIndex) {
        this.index = index;
        this.referenceSequence = referenceSequence;
        this.bins = bins;
        this.linearIndex = linearIndex;
    }

    public int getReferenceSequence() {
        return referenceSequence;
    }

    /**
     * Perform an overlapping query of all bins bounding the given location.
     * @param bin The bin over which to perform an overlapping query.
     * @return The file pointers
     */
    public GATKBAMFileSpan getSpanOverlapping(final Bin bin) {
        if(bin == null)
            return null;

        GATKBin gatkBin = new GATKBin(bin);

        final int binLevel = index.getLevelForBin(bin);
        final int firstLocusInBin = index.getFirstLocusInBin(bin);

        // Add the specified bin to the tree if it exists.
        List<GATKBin> binTree = new ArrayList<GATKBin>();
        if(gatkBin.getBinNumber() < bins.size() && bins.get(gatkBin.getBinNumber()) != null)
            binTree.add(bins.get(gatkBin.getBinNumber()));

        int currentBinLevel = binLevel;
        while(--currentBinLevel >= 0) {
            final int binStart = index.getFirstBinInLevel(currentBinLevel);
            final int binWidth = index.getMaxAddressibleGenomicLocation()/index.getLevelSize(currentBinLevel);
            final int binNumber = firstLocusInBin/binWidth + binStart;
            if(binNumber < bins.size() && bins.get(binNumber) != null)
                binTree.add(bins.get(binNumber));
        }

        List<GATKChunk> chunkList = new ArrayList<GATKChunk>();
        for(GATKBin coveringBin: binTree) {
            for(GATKChunk chunk: coveringBin.getChunkList())
                chunkList.add(chunk.clone());
        }

        final int start = index.getFirstLocusInBin(bin);
        chunkList = optimizeChunkList(chunkList,linearIndex.getMinimumOffset(start));
        return new GATKBAMFileSpan(chunkList.toArray(new GATKChunk[chunkList.size()]));
    }

    private List<GATKChunk> optimizeChunkList(final List<GATKChunk> chunks, final long minimumOffset) {
        GATKChunk lastChunk = null;
        Collections.sort(chunks);
        final List<GATKChunk> result = new ArrayList<GATKChunk>();
        for (final GATKChunk chunk : chunks) {
            if (chunk.getChunkEnd() <= minimumOffset) {
                continue;               // linear index optimization
            }
            if (result.isEmpty()) {
                result.add(chunk);
                lastChunk = chunk;
                continue;
            }
            // Coalesce chunks that are in adjacent file blocks.
            // This is a performance optimization.
            if (!lastChunk.overlaps(chunk) && !lastChunk.isAdjacentTo(chunk)) {
                result.add(chunk);
                lastChunk = chunk;
            } else {
                if (chunk.getChunkEnd() > lastChunk.getChunkEnd()) {
                    lastChunk.setChunkEnd(chunk.getChunkEnd());
                }
            }
        }
        return result;
    }


}
