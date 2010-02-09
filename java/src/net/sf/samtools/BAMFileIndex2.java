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
package net.sf.samtools;


import net.sf.samtools.util.RuntimeIOException;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.util.*;

/**
 * Class for reading BAM file indexes.
 */
public class BAMFileIndex2 extends BAMFileIndex
{
    private static final int MAX_BINS = 37450; // =(8^6-1)/7+1
    private static final int BAM_LIDX_SHIFT = 14;

    /**
     * What is the starting bin for each level?
     */
    private static final int[] LEVEL_STARTS = {0,1,9,73,585,4681};

    /**
     * A mapping of reference sequence index to list of bins.
     */
    protected final SortedMap<Integer,Bin[]> referenceToBins = new TreeMap<Integer,Bin[]>();

    /**
     * A mapping of reference sequence index to linear indices.
     */
    protected final SortedMap<Integer,LinearIndex> referenceToLinearIndices = new TreeMap<Integer,LinearIndex>();

    /**
     * A mapping from bin to the chunks contained in that bin.
     */
    protected final SortedMap<Bin,List<Chunk>> binToChunks = new TreeMap<Bin,List<Chunk>>();

    protected BAMFileIndex2(final File file) {
        super(file);
        loadIndex(file);
    }

    /**
     * Get the number of levels employed by this index.
     * @return Number of levels in this index.
     */
    protected int getNumIndexLevels() {
        return LEVEL_STARTS.length;
    }

    /**
     * Gets the level associated with the given bin number.
     * @param binNumber The bin number for which to determine the level.
     * @return the level associated with the given bin number.
     */
    protected int getLevelForBinNumber(final int binNumber) {
        if(binNumber >= MAX_BINS)
            throw new SAMException("Tried to get level for invalid bin.");
        for(int i = getNumIndexLevels()-1; i >= 0; i--) {
            if(binNumber >= LEVEL_STARTS[i])
                return i;
        }
        throw new SAMException("Unable to find correct bin for bin number "+binNumber);
    }

    /**
     * Completely load the index into memory.
     * @param file File to load.
     */
    private void loadIndex(final File file)  {
        FileInputStream fileStream;
        FileChannel fileChannel;
        MappedByteBuffer fileBuffer;

        try {
            fileStream = new FileInputStream(file);
            fileChannel =  fileStream.getChannel();
            fileBuffer = fileChannel.map(FileChannel.MapMode.READ_ONLY, 0L, fileChannel.size());
            fileBuffer.order(ByteOrder.LITTLE_ENDIAN);
        } catch (IOException exc) {
            throw new RuntimeIOException(exc.getMessage(), exc);
        }

        try {
            final byte[] buffer = new byte[4];
            readBytes(fileBuffer,buffer);
            if (!Arrays.equals(buffer, BAMFileConstants.BAM_INDEX_MAGIC)) {
                throw new RuntimeException("Invalid file header in BAM index " + file +
                        ": " + new String(buffer));
            }

            final int sequenceCount = readInteger(fileBuffer);
            for(int sequence = 0; sequence < sequenceCount; sequence++) {
                final int binCount = readInteger(fileBuffer);
                final Bin[] bins = new Bin[binCount];
                for(int bin = 0; bin < binCount; bin++) {
                    List<Chunk> chunkList = new ArrayList<Chunk>();
                    final int indexBin = readInteger(fileBuffer);
                    final int nChunks = readInteger(fileBuffer);
                    for (int ci = 0; ci < nChunks; ci++) {
                        final long chunkBegin = readLong(fileBuffer);
                        final long chunkEnd = readLong(fileBuffer);
                        chunkList.add(new Chunk(chunkBegin, chunkEnd));
                    }
                    bins[bin] = new Bin(sequence,indexBin);
                    binToChunks.put(bins[bin],chunkList);
                }
                referenceToBins.put(sequence,bins);

                int linearIndexSize = readInteger(fileBuffer);
                long[] linearIndex = new long[linearIndexSize];
                for(int indexEntry = 0; indexEntry < linearIndexSize; indexEntry++)
                    linearIndex[indexEntry] = readLong(fileBuffer);

                referenceToLinearIndices.put(sequence,new LinearIndex(sequence,linearIndex));
            }
        }
        finally {
            try {
                fileChannel.close();
                fileStream.close();
            } catch (IOException exc) {
                throw new RuntimeIOException(exc.getMessage(), exc);
            }
        }
    }

    /**
     * Get list of regions of BAM file that may contain SAMRecords for the given range
     * @param referenceIndex sequence of desired SAMRecords
     * @param startPos 1-based start of the desired interval, inclusive
     * @param endPos 1-based end of the desired interval, inclusive
     * @return array of pairs of virtual file positions.  Each pair is the first and last
     * virtual file position in a range that can be scanned to find SAMRecords that overlap the given
     * positions. The last position in each pair is a virtual file pointer to the first SAMRecord beyond
     * the range that may contain the indicated SAMRecords.
     */
    long[] getFilePointersContaining(final int referenceIndex, final int startPos, final int endPos) {
        List<Bin> bins = getBinsContaining(referenceIndex,startPos,endPos);
        // System.out.println("# Sequence target TID: " + referenceIndex);
        if (bins.size() == 0) {
            return null;
        }

        List<Chunk> chunkList = new ArrayList<Chunk>();
        for(Bin bin: bins)
            chunkList.addAll(binToChunks.get(bin));

        if (chunkList.isEmpty()) {
            return null;
        }

        final int start = (startPos <= 0) ? 0 : startPos-1;
        final int regionLinearBin = start >> BAM_LIDX_SHIFT;
        // System.out.println("# regionLinearBin: " + regionLinearBin);
        LinearIndex index = referenceToLinearIndices.get(referenceIndex);
        long minimumOffset = 0;
        if (regionLinearBin < index.indexEntries.length)
            minimumOffset = index.indexEntries[regionLinearBin];
        chunkList = optimizeChunkList(chunkList, minimumOffset);
        return convertToArray(chunkList);
    }

    long[] getFilePointersBounding(final Bin bin) {
        return convertToArray(binToChunks.get(bin));
    }

    /**
     * Get a list of bins in the BAM file that may contain SAMRecords for the given range.
     * @param referenceIndex sequence of desired SAMRecords
     * @param startPos 1-based start of the desired interval, inclusive
     * @param endPos 1-based end of the desired interval, inclusive
     * @return a list of bins that contain relevant data.
     */
    List<Bin> getBinsContaining(final int referenceIndex, final int startPos, final int endPos) {
        List<Bin> filteredBins = new ArrayList<Bin>();

        if (referenceIndex >= referenceToBins.size()) {
            return null;
        }

        final BitSet regionBins = regionToBins(startPos, endPos);
        if (regionBins == null) {
            return null;
        }

        Bin[] bins = referenceToBins.get(referenceIndex);

        for(Bin bin: bins) {
            if (regionBins.get(bin.binNumber))
                filteredBins.add(bin);
        }

        return filteredBins;
    }

    /**
     * Use to get close to the unmapped reads at the end of a BAM file.
     * @return The file offset of the first record in the last linear bin, or -1
     * if there are no elements in linear bins (i.e. no mapped reads).
     */
    long getStartOfLastLinearBin() {
        LinearIndex lastLinearIndex = referenceToLinearIndices.get(referenceToLinearIndices.lastKey());
        return lastLinearIndex.indexEntries[lastLinearIndex.indexEntries.length-1];
    }

    private List<Chunk> optimizeChunkList(final List<Chunk> chunkList, final long minimumOffset) {
        Chunk lastChunk = null;
        Collections.sort(chunkList);
        final List<Chunk> result = new ArrayList<Chunk>();
        for (final Chunk chunk : chunkList) {
            if (chunk.getChunkEnd() <= minimumOffset) {
                continue;
            }
            if (result.isEmpty()) {
                result.add(chunk);
                lastChunk = chunk;
                continue;
            }
            // Coalesce chunks that are in adjacent file blocks.
            // This is a performance optimization.
            final long lastFileBlock = getFileBlock(lastChunk.getChunkEnd());
            final long chunkFileBlock = getFileBlock(chunk.getChunkStart());
            if (chunkFileBlock - lastFileBlock > 1) {
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

    private long[] convertToArray(final List<Chunk> chunkList) {
        final int count = chunkList.size() * 2;
        if (count == 0) {
            return null;
        }
        int index = 0;
        final long[] result = new long[count];
        for (final Chunk chunk : chunkList) {
            result[index++] = chunk.getChunkStart();
            result[index++] = chunk.getChunkEnd();
        }
        return result;
    }

    /**
     * Get candidate bins for the specified region
     * @param startPos 1-based start of target region, inclusive.
     * @param endPos 1-based end of target region, inclusive.
     * @return bit set for each bin that may contain SAMRecords in the target region.
     */
    protected BitSet regionToBins(final int startPos, final int endPos) {
        final int maxPos = 0x1FFFFFFF;
        final int start = (startPos <= 0) ? 0 : (startPos-1) & maxPos;
        final int end = (endPos <= 0) ? maxPos : (endPos-1) & maxPos;
        if (start > end) {
            return null;
        }
        int k;
        final BitSet bitSet = new BitSet(MAX_BINS);
        bitSet.set(0);
        for (k = LEVEL_STARTS[1] + (start>>26); k <= LEVEL_STARTS[1] + (end>>26); ++k) bitSet.set(k);
        for (k = LEVEL_STARTS[2] + (start>>23); k <= LEVEL_STARTS[2] + (end>>23); ++k) bitSet.set(k);
        for (k = LEVEL_STARTS[3] + (start>>20); k <= LEVEL_STARTS[3] + (end>>20); ++k) bitSet.set(k);
        for (k = LEVEL_STARTS[4] + (start>>17); k <= LEVEL_STARTS[4] + (end>>17); ++k) bitSet.set(k);
        for (k = LEVEL_STARTS[5] + (start>>14); k <= LEVEL_STARTS[5] + (end>>14); ++k) bitSet.set(k);
        return bitSet;
    }

    private long getFileBlock(final long bgzfOffset) {
        return ((bgzfOffset >> 16L) & 0xFFFFFFFFFFFFL);
    }

    private void readBytes(MappedByteBuffer source, final byte[] target) {
        source.get(target);
    }

    private int readInteger(MappedByteBuffer source) {
        return source.getInt();
    }

    private long readLong(MappedByteBuffer source) {
        return source.getLong();
    }
}
