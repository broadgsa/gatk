/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.samtools.Bin;
import net.sf.samtools.GATKBAMFileSpan;
import net.sf.samtools.GATKBin;
import net.sf.samtools.GATKChunk;
import net.sf.samtools.LinearIndex;
import net.sf.samtools.SAMException;
import net.sf.samtools.util.RuntimeIOException;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.*;

/**
 * A basic interface for querying BAM indices.
 *
 * @author mhanna
 * @version 0.1
 */
public class GATKBAMIndex {
    /**
     * Reports the total amount of genomic data that any bin can index.
     */
    protected static final int BIN_GENOMIC_SPAN = 512*1024*1024;

    /**
     * What is the starting bin for each level?
     */
    private static final int[] LEVEL_STARTS = {0,1,9,73,585,4681};

    /**
     * Reports the maximum number of bins that can appear in a BAM file.
     */
    public static final int MAX_BINS = 37450;   // =(8^6-1)/7+1

    private final File mFile;

    private final List<GATKBin> bins = new ArrayList<GATKBin>();
    private final LinearIndex linearIndex;

    public GATKBAMIndex(final File file, final int referenceSequence) {
        mFile = file;
        // Open the file stream.

        FileInputStream fileStream = null;
        FileChannel fileChannel  = null;

        try {
            fileStream = new FileInputStream(mFile);
            fileChannel = fileStream.getChannel();
        }
        catch (IOException exc) {
            throw new RuntimeIOException(exc.getMessage(), exc);
        }

        // Verify the magic number.
        seek(fileChannel,0);
        final byte[] buffer = new byte[4];
        readBytes(fileChannel,buffer);
        if (!Arrays.equals(buffer, GATKBAMFileConstants.BAM_INDEX_MAGIC)) {
            throw new RuntimeException("Invalid file header in BAM index " + mFile +
                                       ": " + new String(buffer));
        }

        seek(fileChannel,4);

        final int sequenceCount = readInteger(fileChannel);

        if (referenceSequence >= sequenceCount)
            throw new ReviewedStingException("Invalid sequence number " + referenceSequence);

        skipToSequence(fileChannel,referenceSequence);

        int binCount = readInteger(fileChannel);
        for (int binNumber = 0; binNumber < binCount; binNumber++) {
            final int indexBin = readInteger(fileChannel);
            final int nChunks = readInteger(fileChannel);
            List<GATKChunk> chunks = new ArrayList<GATKChunk>(nChunks);
            for (int ci = 0; ci < nChunks; ci++) {
                final long chunkBegin = readLong(fileChannel);
                final long chunkEnd = readLong(fileChannel);
                chunks.add(new GATKChunk(chunkBegin, chunkEnd));
            }
            GATKBin bin = new GATKBin(referenceSequence, indexBin);
            bin.setChunkList(chunks.toArray(new GATKChunk[chunks.size()]));
            while(indexBin >= bins.size())
                bins.add(null);
            bins.set(indexBin,bin);
        }

        final int nLinearBins = readInteger(fileChannel);
        long[] linearIndexEntries = new long[nLinearBins];
        for(int linearIndexOffset = 0; linearIndexOffset < nLinearBins; linearIndexOffset++)
            linearIndexEntries[linearIndexOffset] = readLong(fileChannel);

        linearIndex = new LinearIndex(referenceSequence,0,linearIndexEntries);

        try {
            fileChannel.close();
            fileStream.close();
        }
        catch (IOException exc) {
            throw new RuntimeIOException(exc.getMessage(), exc);
        }
    }

    /**
     * Get the number of levels employed by this index.
     * @return Number of levels in this index.
     */
    public static int getNumIndexLevels() {
        return LEVEL_STARTS.length;
    }

    /**
     * Gets the first bin in the given level.
     * @param levelNumber Level number.  0-based.
     * @return The first bin in this level.
     */
    public static int getFirstBinInLevel(final int levelNumber) {
        return LEVEL_STARTS[levelNumber];
    }

    /**
     * Gets the number of bins in the given level.
     * @param levelNumber Level number.  0-based.
     * @return The size (number of possible bins) of the given level.
     */
    public int getLevelSize(final int levelNumber) {
        if(levelNumber == getNumIndexLevels()-1)
            return MAX_BINS-LEVEL_STARTS[levelNumber]-1;
        else
            return LEVEL_STARTS[levelNumber+1]-LEVEL_STARTS[levelNumber];
    }

    /**
     * Gets the level associated with the given bin number.
     * @param bin The bin  for which to determine the level.
     * @return the level associated with the given bin number.
     */
    public int getLevelForBin(final Bin bin) {
        GATKBin gatkBin = new GATKBin(bin);
        if(gatkBin.getBinNumber() >= MAX_BINS)
            throw new SAMException("Tried to get level for invalid bin.");
        for(int i = getNumIndexLevels()-1; i >= 0; i--) {
            if(gatkBin.getBinNumber() >= LEVEL_STARTS[i])
                return i;
        }
        throw new SAMException("Unable to find correct bin for bin "+bin);
    }

    /**
     * Gets the first locus that this bin can index into.
     * @param bin The bin to test.
     * @return The last position that the given bin can represent.
     */
    public int getFirstLocusInBin(final Bin bin) {
        final int level = getLevelForBin(bin);
        final int levelStart = LEVEL_STARTS[level];
        final int levelSize = ((level==getNumIndexLevels()-1) ? MAX_BINS-1 : LEVEL_STARTS[level+1]) - levelStart;
        return (new GATKBin(bin).getBinNumber() - levelStart)*(BIN_GENOMIC_SPAN /levelSize)+1;
    }

    /**
     * Gets the last locus that this bin can index into.
     * @param bin The bin to test.
     * @return The last position that the given bin can represent.
     */
    public int getLastLocusInBin(final Bin bin) {
        final int level = getLevelForBin(bin);
        final int levelStart = LEVEL_STARTS[level];
        final int levelSize = ((level==getNumIndexLevels()-1) ? MAX_BINS-1 : LEVEL_STARTS[level+1]) - levelStart;
        return (new GATKBin(bin).getBinNumber()-levelStart+1)*(BIN_GENOMIC_SPAN /levelSize);
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

        final int binLevel = getLevelForBin(bin);
        final int firstLocusInBin = getFirstLocusInBin(bin);

        // Add the specified bin to the tree if it exists.
        List<GATKBin> binTree = new ArrayList<GATKBin>();
        if(gatkBin.getBinNumber() < bins.size() && bins.get(gatkBin.getBinNumber()) != null)
            binTree.add(bins.get(gatkBin.getBinNumber()));

        int currentBinLevel = binLevel;
        while(--currentBinLevel >= 0) {
            final int binStart = getFirstBinInLevel(currentBinLevel);
            final int binWidth = getMaxAddressibleGenomicLocation()/getLevelSize(currentBinLevel);
            final int binNumber = firstLocusInBin/binWidth + binStart;
            if(binNumber < bins.size() && bins.get(binNumber) != null)
                binTree.add(bins.get(binNumber));
        }

        List<GATKChunk> chunkList = new ArrayList<GATKChunk>();
        for(GATKBin coveringBin: binTree) {
            for(GATKChunk chunk: coveringBin.getChunkList())
                chunkList.add(chunk.clone());
        }

        final int start = getFirstLocusInBin(bin);
        chunkList = optimizeChunkList(chunkList,linearIndex.getMinimumOffset(start));
        return new GATKBAMFileSpan(chunkList.toArray(new GATKChunk[chunkList.size()]));
    }

    protected List<GATKChunk> optimizeChunkList(final List<GATKChunk> chunks, final long minimumOffset) {
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

    /**
     * Gets the possible number of bins for a given reference sequence.
     * @return How many bins could possibly be used according to this indexing scheme to index a single contig.
     */
    protected int getMaxAddressibleGenomicLocation() {
        return BIN_GENOMIC_SPAN;
    }    

    protected void skipToSequence(final FileChannel fileChannel, final int sequenceIndex) {
        for (int i = 0; i < sequenceIndex; i++) {
            // System.out.println("# Sequence TID: " + i);
            final int nBins = readInteger(fileChannel);
            // System.out.println("# nBins: " + nBins);
            for (int j = 0; j < nBins; j++) {
                final int bin = readInteger(fileChannel);
                final int nChunks = readInteger(fileChannel);
                // System.out.println("# bin[" + j + "] = " + bin + ", nChunks = " + nChunks);
                skipBytes(fileChannel, 16 * nChunks);
            }
            final int nLinearBins = readInteger(fileChannel);
            // System.out.println("# nLinearBins: " + nLinearBins);
            skipBytes(fileChannel, 8 * nLinearBins);
        }
    }

    private void readBytes(final FileChannel fileChannel, final byte[] bytes) {
        ByteBuffer buffer = ByteBuffer.wrap(bytes);
        try {
            fileChannel.read(buffer);
        }
        catch(IOException ex) {
            throw new ReviewedStingException("Index: unable to read bytes from index file.");
        }
    }

    private static final int INT_SIZE_IN_BYTES = Integer.SIZE / 8;
    private static final int LONG_SIZE_IN_BYTES = Long.SIZE / 8;

    private ByteBuffer wrapBuffer(final byte[] bytes) {
        ByteBuffer buffer = ByteBuffer.wrap(bytes);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        return buffer;
    }

    private int readInteger(final FileChannel fileChannel) {
        byte[] bytes = new byte[INT_SIZE_IN_BYTES];
        readBytes(fileChannel,bytes);
        return wrapBuffer(bytes).getInt();
    }

    private long readLong(final FileChannel fileChannel) {
        byte[] bytes = new byte[LONG_SIZE_IN_BYTES];
        readBytes(fileChannel,bytes);
        return wrapBuffer(bytes).getLong();
    }

    private void skipBytes(final FileChannel fileChannel, final int count) {
        try {
            fileChannel.position(fileChannel.position() + count);
        }
        catch(IOException ex) {
            throw new ReviewedStingException("Index: unable to reposition file channel.");
        }
    }

    private void seek(final FileChannel fileChannel, final int position) {
        try {
            fileChannel.position(position);
        }
        catch(IOException ex) {
            throw new ReviewedStingException("Index: unable to reposition of file channel.");
        }
    }
}
