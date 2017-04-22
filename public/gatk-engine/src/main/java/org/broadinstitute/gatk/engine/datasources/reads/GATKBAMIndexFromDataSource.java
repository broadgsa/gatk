/*
* Copyright 2012-2016 Broad Institute, Inc.
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
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * This class implements BAM index querying API
 * by wrapping a class BrowseableBAMIndex from HTSJDK.
 *
 * @version 0.1
 */
public class GATKBAMIndexFromDataSource extends GATKBAMIndex {
    private File sourceFile;
    private SAMFileHeader sourceHeader;
    private BrowseableBAMIndex index;

    public GATKBAMIndexFromDataSource(final File sourceFile, final SAMFileHeader sourceHeader, final BrowseableBAMIndex index) {
        this.sourceFile = sourceFile;
        this.sourceHeader = sourceHeader;
        this.index = index;
    }

    @Override
    public GATKBAMIndexData readReferenceSequence(final int referenceSequence) {
        final List<SAMSequenceRecord> sequences = sourceHeader.getSequenceDictionary().getSequences();
        if (referenceSequence >= sequences.size())
            throw new ReviewedGATKException("Sequence number " + referenceSequence + " cannot be greater or equal to " + sequences.size() + " in index file " + sourceFile);


        final BinList sourceBins = index.getBinsOverlapping(referenceSequence, 0, sequences.get(referenceSequence).getSequenceLength());

        final List<GATKBin> bins = new ArrayList<>();
        for (Bin sourceBin : sourceBins) {
            final int indexBin = sourceBin.getBinNumber();
            while(indexBin >= bins.size())
                bins.add(null);

            final GATKBin bin = new GATKBin(referenceSequence, indexBin);
            final List<Chunk> chunks = index.getSpanOverlapping(sourceBin).getChunks();
            final List<GATKChunk> gatkChunks = new ArrayList<>(chunks.size());
            for (Chunk chunk : chunks) {
                gatkChunks.add(new GATKChunk(chunk));
            }

            bin.setChunkList(gatkChunks.toArray(new GATKChunk[gatkChunks.size()]));

            bins.set(indexBin, bin);
        }

        // there is no interface to get linear index from HTSJDK
        final LinearIndex linearIndex = new LinearIndex(referenceSequence, 0, new long[]{});

        return new GATKBAMIndexData(this,referenceSequence,bins,linearIndex);
    }

    @Override
    public int getLevelSize(int levelNumber) {
        return index.getLevelSize(levelNumber);
    }

    @Override
    public int getLevelForBin(Bin bin) {
        return index.getLevelForBin(bin);
    }

    @Override
    public int getFirstLocusInBin(Bin bin) {
        return index.getFirstLocusInBin(bin);
    }

    @Override
    public int getLastLocusInBin(Bin bin) {
        return index.getLastLocusInBin(bin);
    }

    @Override
    public long getStartOfLastLinearBin() {
        return index.getStartOfLastLinearBin();
    }
}
