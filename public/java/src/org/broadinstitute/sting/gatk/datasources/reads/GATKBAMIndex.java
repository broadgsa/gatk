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

import net.sf.samtools.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A basic interface for querying BAM indices.
 * Very much not thread-safe.
 *
 * @author mhanna
 * @version 0.1
 */
public class GATKBAMIndex {
    /**
     * BAM index file magic number.
     */
    private static final byte[] BAM_INDEX_MAGIC = "BAI\1".getBytes();

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

    /**
     * Number of sequences stored in this index.
     */
    private final int sequenceCount;

    /**
     * A cache of the starting positions of the sequences.
     */
    private final long[] sequenceStartCache;

    private FileInputStream fileStream;
    private FileChannel fileChannel;

    public GATKBAMIndex(final File file) {
        mFile = file;
        // Open the file stream.
        openIndexFile();

        // Verify the magic number.
        seek(0);
        final byte[] buffer = readBytes(4);
        if (!Arrays.equals(buffer, BAM_INDEX_MAGIC)) {
            throw new RuntimeException("Invalid file header in BAM index " + mFile +
                                       ": " + new String(buffer));
        }

        seek(4);

        sequenceCount = readInteger();

        // Create a cache of the starting position of each sequence.  Initialize it to -1.
        sequenceStartCache = new long[sequenceCount];
        for(int i = 1; i < sequenceCount; i++)
            sequenceStartCache[i] = -1;

        // Seed the first element in the array with the current position.
        if(sequenceCount > 0)
            sequenceStartCache[0] = position();

        closeIndexFile();
    }

    public GATKBAMIndexData readReferenceSequence(final int referenceSequence) {
        openIndexFile();

        if (referenceSequence >= sequenceCount)
            throw new ReviewedStingException("Invalid sequence number " + referenceSequence);

        skipToSequence(referenceSequence);

        int binCount = readInteger();
        List<GATKBin> bins = new ArrayList<GATKBin>();
        for (int binNumber = 0; binNumber < binCount; binNumber++) {
            final int indexBin = readInteger();
            final int nChunks = readInteger();

            List<GATKChunk> chunks = new ArrayList<GATKChunk>(nChunks);
            long[] rawChunkData = readLongs(nChunks*2);
            for (int ci = 0; ci < nChunks; ci++) {
                final long chunkBegin = rawChunkData[ci*2];
                final long chunkEnd = rawChunkData[ci*2+1];
                chunks.add(new GATKChunk(chunkBegin, chunkEnd));
            }
            GATKBin bin = new GATKBin(referenceSequence, indexBin);
            bin.setChunkList(chunks.toArray(new GATKChunk[chunks.size()]));
            while(indexBin >= bins.size())
                bins.add(null);
            bins.set(indexBin,bin);
        }

        final int nLinearBins = readInteger();
        long[] linearIndexEntries = readLongs(nLinearBins);

        LinearIndex linearIndex = new LinearIndex(referenceSequence,0,linearIndexEntries);

        closeIndexFile();

        return new GATKBAMIndexData(this,referenceSequence,bins,linearIndex);
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
     * Use to get close to the unmapped reads at the end of a BAM file.
     * @return The file offset of the first record in the last linear bin, or -1
     * if there are no elements in linear bins (i.e. no mapped reads).
     */
    public long getStartOfLastLinearBin() {
        openIndexFile();

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
                lastLinearIndexPointer = readLongs(1)[0];
            }
        }

        closeIndexFile();

        return lastLinearIndexPointer;
    }

    /**
     * Gets the possible number of bins for a given reference sequence.
     * @return How many bins could possibly be used according to this indexing scheme to index a single contig.
     */
    protected int getMaxAddressibleGenomicLocation() {
        return BIN_GENOMIC_SPAN;
    }    

    protected void skipToSequence(final int referenceSequence) {
        // Find the offset in the file of the last sequence whose position has been determined.  Start here
        // when searching the sequence for the next value to read.  (Note that sequenceStartCache[0] will always
        // be present, so no extra stopping condition is necessary.
        int sequenceIndex = referenceSequence;
        while(sequenceStartCache[sequenceIndex] == -1)
            sequenceIndex--;

        // Advance to the most recently found position.
        seek(sequenceStartCache[sequenceIndex]);

        for (int i = sequenceIndex; i < referenceSequence; i++) {
            sequenceStartCache[i] = position();

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

        sequenceStartCache[referenceSequence] = position();
    }

    private void openIndexFile() {
        try {
            fileStream = new FileInputStream(mFile);
            fileChannel = fileStream.getChannel();
        }
        catch (IOException exc) {
            throw new ReviewedStingException("Unable to open index file " + mFile, exc);            
        }
    }

    private void closeIndexFile() {
        try {
            fileChannel.close();
            fileStream.close();
        }
        catch (IOException exc) {
            throw new ReviewedStingException("Unable to close index file " + mFile, exc);
        }
    }

    private static final int INT_SIZE_IN_BYTES = Integer.SIZE / 8;
    private static final int LONG_SIZE_IN_BYTES = Long.SIZE / 8;

    private byte[] readBytes(int count) {
        ByteBuffer buffer = getBuffer(count);
        read(buffer);
        buffer.flip();
        byte[] contents = new byte[count];
        buffer.get(contents);
        return contents;
    }

    private int readInteger() {
        ByteBuffer buffer = getBuffer(INT_SIZE_IN_BYTES);
        read(buffer);
        buffer.flip();
        return buffer.getInt();
    }

    /**
     * Reads an array of <count> longs from the file channel, returning the results as an array.
     * @param count Number of longs to read.
     * @return An array of longs.  Size of array should match count.
     */
    private long[] readLongs(final int count) {
        ByteBuffer buffer = getBuffer(count*LONG_SIZE_IN_BYTES);
        read(buffer);
        buffer.flip();
        long[] result = new long[count];
        for(int i = 0; i < count; i++)
            result[i] = buffer.getLong();
        return result;
    }

    private void read(final ByteBuffer buffer) {
        try {
            fileChannel.read(buffer);
        }
        catch(IOException ex) {
            throw new ReviewedStingException("Index: unable to read bytes from index file.");
        }
    }


    /**
     * A reusable buffer for use by this index generator.
     * TODO: Should this be a SoftReference?
     */
    private ByteBuffer buffer = null;

    private ByteBuffer getBuffer(final int size) {
        if(buffer == null || buffer.capacity() < size) {
            // Allocate a new byte buffer.  For now, make it indirect to make sure it winds up on the heap for easier debugging.
            buffer = ByteBuffer.allocate(size);
            buffer.order(ByteOrder.LITTLE_ENDIAN);
        }
        buffer.clear();
        buffer.limit(size);
        return buffer;
    }

    private void skipBytes(final int count) {
        try {
            fileChannel.position(fileChannel.position() + count);
        }
        catch(IOException ex) {
            throw new ReviewedStingException("Index: unable to reposition file channel.");
        }
    }

    private void seek(final long position) {
        try {
            fileChannel.position(position);
        }
        catch(IOException ex) {
            throw new ReviewedStingException("Index: unable to reposition of file channel.");
        }
    }

    /**
     * Retrieve the position from the current file channel.
     * @return position of the current file channel.
     */
    private long position() {
        try {
            return fileChannel.position();
        }
        catch (IOException exc) {
            throw new ReviewedStingException("Unable to read position from index file " + mFile, exc);
        }
    }    
}
