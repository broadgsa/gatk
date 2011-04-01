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

import net.sf.samtools.BAMIndex;
import net.sf.samtools.BAMIndexMetaData;
import net.sf.samtools.Bin;
import net.sf.samtools.BrowseableBAMIndex;
import net.sf.samtools.GATKBAMFileSpan;
import net.sf.samtools.GATKBin;
import net.sf.samtools.GATKBinList;
import net.sf.samtools.GATKChunk;
import net.sf.samtools.LinearIndex;
import net.sf.samtools.SAMException;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.*;

/**
 * A basic interface for querying BAM indices.
 *
 * @author mhanna
 * @version 0.1
 */
public class GATKBAMIndex implements BAMIndex, BrowseableBAMIndex {
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

    public static final int MAX_LINEAR_INDEX_SIZE = MAX_BINS+1-LEVEL_STARTS[LEVEL_STARTS.length-1];

    private final File mFile;
    private final MappedByteBuffer mFileBuffer;

    private SAMSequenceDictionary mBamDictionary = null;

    public GATKBAMIndex(final File file, final SAMSequenceDictionary dictionary) {
        mFile = file;
        mBamDictionary = dictionary;
        // Open the file stream.
        try {
            FileInputStream fileStream = new FileInputStream(mFile);
            FileChannel fileChannel = fileStream.getChannel();
            mFileBuffer = fileChannel.map(FileChannel.MapMode.READ_ONLY, 0L, fileChannel.size());
            mFileBuffer.order(ByteOrder.LITTLE_ENDIAN);

            fileChannel.close();
            fileStream.close();
        }
        catch (IOException exc) {
            throw new RuntimeIOException(exc.getMessage(), exc);
        }

        // Verify the magic number.
        seek(0);
        final byte[] buffer = new byte[4];
        readBytes(buffer);
        if (!Arrays.equals(buffer, GATKBAMFileConstants.BAM_INDEX_MAGIC)) {
            throw new RuntimeException("Invalid file header in BAM index " + mFile +
                                       ": " + new String(buffer));
        }
    }

    /**
     * Gets the file backing this index.
     * @return The index file.
     */
    public File getIndexFile() {
        return mFile;
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
    @Override
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
    @Override
    public int getLastLocusInBin(final Bin bin) {
        final int level = getLevelForBin(bin);
        final int levelStart = LEVEL_STARTS[level];
        final int levelSize = ((level==getNumIndexLevels()-1) ? MAX_BINS-1 : LEVEL_STARTS[level+1]) - levelStart;
        return (new GATKBin(bin).getBinNumber()-levelStart+1)*(BIN_GENOMIC_SPAN /levelSize);
    }

    public int getNumberOfReferences() {
        seek(4);
        return readInteger();
    }

    /**
     * Use to get close to the unmapped reads at the end of a BAM file.
     * @return The file offset of the first record in the last linear bin, or -1
     * if there are no elements in linear bins (i.e. no mapped reads).
     */
    public long getStartOfLastLinearBin() {
        seek(4);

        final int sequenceCount = readInteger();
        // Because no reads may align to the last sequence in the sequence dictionary,
        // grab the last element of the linear index for each sequence, and return
        // the last one from the last sequence that has one.
        long lastLinearIndexPointer = -1;
        for (int i = 0; i < sequenceCount; i++) {
            // System.out.println("# Sequence TID: " + i);
            final int nBins = readInteger();
            // System.out.println("# nBins: " + nBins);
            for (int j1 = 0; j1 < nBins; j1++) {
                // Skip bin #
                skipBytes(4);
                final int nChunks = readInteger();
                // Skip chunks
                skipBytes(16 * nChunks);
            }
            final int nLinearBins = readInteger();
            if (nLinearBins > 0) {
                // Skip to last element of list of linear bins
                skipBytes(8 * (nLinearBins - 1));
                lastLinearIndexPointer = readLong();
            }
        }

        return lastLinearIndexPointer;
    }

    /**
     * Gets meta data for the given reference including information about number of aligned, unaligned, and noCoordinate records
     * @param reference the reference of interest
     * @return meta data for the reference
     */
    public BAMIndexMetaData getMetaData(int reference) {
        throw new UnsupportedOperationException("Cannot retrieve metadata for GATKBAMIndex");
    }

    /**
     * Returns count of records unassociated with any reference. Call before the index file is closed
     *
     * @return meta data at the end of the bam index that indicates count of records holding no coordinates
     * or null if no meta data (old index format)
     */
    public Long getNoCoordinateCount() {

        seek(4);
        final int sequenceCount = readInteger();

        skipToSequence(sequenceCount);
        try { // in case of old index file without meta data
            return readLong();
        } catch (Exception e) {
            return null;
        }
    }

    /**
     * Get list of regions of BAM file that may contain SAMRecords for the given range
     * @param referenceIndex sequence of desired SAMRecords
     * @param startPos 1-based start of the desired interval, inclusive
     * @param endPos 1-based end of the desired interval, inclusive
     * @return the virtual file position.  Each pair is the first and last virtual file position
     *         in a range that can be scanned to find SAMRecords that overlap the given positions.
     */
    @Override
    public GATKBAMFileSpan getSpanOverlapping(final int referenceIndex, final int startPos, final int endPos) {
        BAMIndexContent queryResults = getQueryResults(referenceIndex);

        if(queryResults == null)
            return null;

        GATKBinList overlappingBins = getBinsOverlapping(referenceIndex,startPos,endPos);

        // System.out.println("# Sequence target TID: " + referenceIndex);
        List<GATKBin> bins = new ArrayList<GATKBin>();
        for(GATKBin bin: queryResults.getBins()) {
            if (overlappingBins.getBins().get(bin.getBinNumber()))
                bins.add(bin);
        }

        if (bins.isEmpty()) {
            return null;
        }

        List<GATKChunk> chunkList = new ArrayList<GATKChunk>();
        for(GATKBin bin: bins) {
            for(GATKChunk chunk: bin.getChunkList())
                chunkList.add(chunk.clone());
        }

        if (chunkList.isEmpty()) {
            return null;
        }

        chunkList = optimizeChunkList(chunkList,queryResults.getLinearIndex().getMinimumOffset(startPos));
        return new GATKBAMFileSpan(chunkList.toArray(new GATKChunk[chunkList.size()]));
    }

    /**
     * Perform an overlapping query of all bins bounding the given location.
     * @param bin The bin over which to perform an overlapping query.
     * @return The file pointers
     */
    @Override
    public GATKBAMFileSpan getSpanOverlapping(final Bin bin) {
        if(bin == null)
            return null;

        GATKBin gatkBin = new GATKBin(bin);

        final int referenceSequence = gatkBin.getReferenceSequence();
        BAMIndexContent indexQuery = getQueryResults(referenceSequence);

        if(indexQuery == null)
            return null;

        final int binLevel = getLevelForBin(bin);
        final int firstLocusInBin = getFirstLocusInBin(bin);

        // Add the specified bin to the tree if it exists.
        List<GATKBin> binTree = new ArrayList<GATKBin>();
        if(indexQuery.containsBin(gatkBin))
            binTree.add(indexQuery.getBins().getBin(gatkBin.getBinNumber()));

        int currentBinLevel = binLevel;
        while(--currentBinLevel >= 0) {
            final int binStart = getFirstBinInLevel(currentBinLevel);
            final int binWidth = getMaxAddressibleGenomicLocation()/getLevelSize(currentBinLevel);
            final int binNumber = firstLocusInBin/binWidth + binStart;
            GATKBin parentBin = indexQuery.getBins().getBin(binNumber);
            if(parentBin != null && indexQuery.containsBin(parentBin))
                binTree.add(parentBin);
        }

        List<GATKChunk> chunkList = new ArrayList<GATKChunk>();
        for(GATKBin coveringBin: binTree) {
            for(GATKChunk chunk: coveringBin.getChunkList())
                chunkList.add(chunk.clone());
        }

        final int start = getFirstLocusInBin(bin);
        chunkList = optimizeChunkList(chunkList,indexQuery.getLinearIndex().getMinimumOffset(start));
        return new GATKBAMFileSpan(chunkList.toArray(new GATKChunk[chunkList.size()]));
    }

    public GATKBAMFileSpan getContentsOfBin(final Bin bin) {
        if(bin == null)
            return null;

        GATKBin gatkBin = new GATKBin(bin);

        BAMIndexContent indexQuery = getQueryResults(gatkBin.getReferenceSequence());

        if(indexQuery == null)
            return null;

        GATKBin queriedBin = indexQuery.getBins().getBin(gatkBin.getBinNumber());

        return queriedBin != null ? new GATKBAMFileSpan(queriedBin.getChunkList()) : null;
    }

    /**
     * Retrieves the linear index for the given reference sequence.
     * @param referenceSequence Reference sequence number for which to retrieve the reference.
     * @return The linear index for the given reference sequence.
     */
    public LinearIndex getLinearIndex(int referenceSequence) {
        return getQueryResults(referenceSequence).getLinearIndex();
    }

    /**
     * Get a list of bins in the BAM file that may contain SAMRecords for the given range.
     * @param referenceIndex sequence of desired SAMRecords
     * @param startPos 1-based start of the desired interval, inclusive
     * @param endPos 1-based end of the desired interval, inclusive
     * @return a list of bins that contain relevant data.
     */
    public GATKBinList getBinsOverlapping(final int referenceIndex, final int startPos, final int endPos) {
        final BitSet regionBins = regionToBins(startPos,endPos);
        if (regionBins == null) {
            return null;
        }
        return new GATKBinList(referenceIndex,regionBins);
    }

    protected BAMIndexContent query(final int referenceSequence, final int startPos, final int endPos) {
        seek(4);

        List<GATKChunk> metaDataChunks = new ArrayList<GATKChunk>();

        final int sequenceCount = readInteger();

        if (referenceSequence >= sequenceCount) {
            return null;
        }

        final BitSet regionBins = regionToBins(startPos, endPos);
        if (regionBins == null) {
            return null;
        }

        skipToSequence(referenceSequence);

        int binCount = readInteger();
        boolean metaDataSeen = false;
        GATKBin[] bins = new GATKBin[getMaxBinNumberForReference(referenceSequence) +1];
        for (int binNumber = 0; binNumber < binCount; binNumber++) {
            final int indexBin = readInteger();
            final int nChunks = readInteger();
            List<GATKChunk> chunks = new ArrayList<GATKChunk>(nChunks);
            // System.out.println("# bin[" + i + "] = " + indexBin + ", nChunks = " + nChunks);
            GATKChunk lastChunk = null;
            if (regionBins.get(indexBin)) {
                for (int ci = 0; ci < nChunks; ci++) {
                    final long chunkBegin = readLong();
                    final long chunkEnd = readLong();
                    lastChunk = new GATKChunk(chunkBegin, chunkEnd);
                    chunks.add(lastChunk);
                }
            } else if (indexBin == MAX_BINS) {
                // meta data - build the bin so that the count of bins is correct;
                // but don't attach meta chunks to the bin, or normal queries will be off
                for (int ci = 0; ci < nChunks; ci++) {
                    final long chunkBegin = readLong();
                    final long chunkEnd = readLong();
                    lastChunk = new GATKChunk(chunkBegin, chunkEnd);
                    metaDataChunks.add(lastChunk);
                }
                metaDataSeen = true;
                continue; // don't create a Bin
            } else {
                skipBytes(16 * nChunks);
            }
            GATKBin bin = new GATKBin(referenceSequence, indexBin);
            bin.setChunkList(chunks.toArray(new GATKChunk[chunks.size()]));
            bins[indexBin] = bin;
        }

        final int nLinearBins = readInteger();

        final int regionLinearBinStart = LinearIndex.convertToLinearIndexOffset(startPos);
        final int regionLinearBinStop = endPos > 0 ? LinearIndex.convertToLinearIndexOffset(endPos) : nLinearBins-1;
        final int actualStop = Math.min(regionLinearBinStop, nLinearBins -1);

        long[] linearIndexEntries = new long[0];
        if (regionLinearBinStart < nLinearBins) {
            linearIndexEntries = new long[actualStop-regionLinearBinStart+1];
            skipBytes(8 * regionLinearBinStart);
            for(int linearBin = regionLinearBinStart; linearBin <= actualStop; linearBin++)
                linearIndexEntries[linearBin-regionLinearBinStart] = readLong();
        }

        final LinearIndex linearIndex = new LinearIndex(referenceSequence,regionLinearBinStart,linearIndexEntries);

        return new BAMIndexContent(referenceSequence, bins, binCount - (metaDataSeen? 1 : 0), linearIndex);
    }

    /**
     * The maxiumum bin number for a reference sequence of a given length
     */
    static int getMaxBinNumberForSequenceLength(int sequenceLength) {
        return getFirstBinInLevel(getNumIndexLevels() - 1) + (sequenceLength >> 14);
        // return 4680 + (sequenceLength >> 14); // note 4680 = getFirstBinInLevel(getNumIndexLevels() - 1)
    }

    /**
     * Looks up the cached BAM query results if they're still in the cache and not expired.  Otherwise,
     * retrieves the cache results from disk.
     * @param referenceIndex The reference to load.  CachingBAMFileIndex only stores index data for entire references.
     * @return The index information for this reference.
     */
    protected BAMIndexContent getQueryResults(final int referenceIndex) {
        // If not in the cache, attempt to load it from disk.
        return query(referenceIndex,1,-1);
    }

    /**
     * Gets the possible number of bins for a given reference sequence.
     * @return How many bins could possibly be used according to this indexing scheme to index a single contig.
     */
    protected int getMaxAddressibleGenomicLocation() {
        return BIN_GENOMIC_SPAN;
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
        for (k =    1 + (start>>26); k <=    1 + (end>>26); ++k) bitSet.set(k);
        for (k =    9 + (start>>23); k <=    9 + (end>>23); ++k) bitSet.set(k);
        for (k =   73 + (start>>20); k <=   73 + (end>>20); ++k) bitSet.set(k);
        for (k =  585 + (start>>17); k <=  585 + (end>>17); ++k) bitSet.set(k);
        for (k = 4681 + (start>>14); k <= 4681 + (end>>14); ++k) bitSet.set(k);
        return bitSet;
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
     * The maximum possible bin number for this reference sequence.
     * This is based on the maximum coordinate position of the reference
     * which is based on the size of the reference
     */
    private int getMaxBinNumberForReference(final int reference) {
        try {
            final int sequenceLength = mBamDictionary.getSequence(reference).getSequenceLength();
            return getMaxBinNumberForSequenceLength(sequenceLength);
        } catch (Exception e) {
            return MAX_BINS;
        }
    }

    protected void skipToSequence(final int sequenceIndex) {
        for (int i = 0; i < sequenceIndex; i++) {
            // System.out.println("# Sequence TID: " + i);
            final int nBins = readInteger();
            // System.out.println("# nBins: " + nBins);
            for (int j = 0; j < nBins; j++) {
                final int bin = readInteger();
                final int nChunks = readInteger();
                // System.out.println("# bin[" + j + "] = " + bin + ", nChunks = " + nChunks);
                skipBytes(16 * nChunks);
            }
            final int nLinearBins = readInteger();
            // System.out.println("# nLinearBins: " + nLinearBins);
            skipBytes(8 * nLinearBins);
        }
    }

    private void readBytes(final byte[] bytes) {
        mFileBuffer.get(bytes);
    }

    protected int readInteger() {
        return mFileBuffer.getInt();
    }

    protected long readLong() {
        return mFileBuffer.getLong();
    }

    protected void skipBytes(final int count) {
        mFileBuffer.position(mFileBuffer.position() + count);
    }

    protected void seek(final int position) {
        mFileBuffer.position(position);
    }

    protected long position() {
        return mFileBuffer.position();
    }

    @Override
    public void close() {
    }
}
