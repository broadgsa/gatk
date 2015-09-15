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
import htsjdk.samtools.seekablestream.SeekableBufferedStream;
import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
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

    private final SAMSequenceDictionary sequenceDictionary;
    private final File mFile;

    //TODO: figure out a good value for this buffer size
    private static final int BUFFERED_STREAM_BUFFER_SIZE = 8192;

    /**
     * Number of sequences stored in this index.
     */
    private final int sequenceCount;

    /**
     * A cache of the starting positions of the sequences.
     */
    private final long[] sequenceStartCache;

    private SeekableFileStream fileStream;
    private SeekableStream baiStream;
    private SeekableBufferedStream bufferedStream;
    private long fileLength;

    public GATKBAMIndex(final File file, final SAMSequenceDictionary sequenceDictionary) {
        mFile = file;
        this.sequenceDictionary = sequenceDictionary;

        // Open the file stream.
        openIndexFile();

        // Verify the magic number.
        seek(0);
        final byte[] buffer = readBytes(4);
        if (!Arrays.equals(buffer, BAM_INDEX_MAGIC)) {
            throw new ReviewedGATKException("Invalid file header in BAM index " + mFile +
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
            throw new ReviewedGATKException("Invalid sequence number " + referenceSequence + " in index file " + mFile);

        skipToSequence(referenceSequence);

        int binCount = readInteger();
        List<GATKBin> bins = new ArrayList<>();
        for (int binNumber = 0; binNumber < binCount; binNumber++) {
            final int indexBin = readInteger();
            final int nChunks = readInteger();

            List<GATKChunk> chunks = new ArrayList<>(nChunks);
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
            throw new ReviewedGATKException("Tried to get level for invalid bin in index file " + mFile);
        for(int i = getNumIndexLevels()-1; i >= 0; i--) {
            if(gatkBin.getBinNumber() >= LEVEL_STARTS[i])
                return i;
        }
        throw new ReviewedGATKException("Unable to find correct bin for bin " + bin + " in index file " + mFile);
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
                /* final int bin = */
                readInteger();
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
            fileStream = new SeekableFileStream(mFile);
            baiStream = SamIndexes.asBaiSeekableStreamOrNull(fileStream, sequenceDictionary);
            bufferedStream = new SeekableBufferedStream(baiStream, BUFFERED_STREAM_BUFFER_SIZE);
            fileLength=bufferedStream.length();
        }
        catch (IOException exc) {
            throw new ReviewedGATKException("Unable to open index file (" + exc.getMessage() +")" + mFile, exc);
        }
    }

    private void closeIndexFile() {
        try {
            bufferedStream.close();
            baiStream.close();
            fileStream.close();
            fileLength = -1;
        }
        catch (IOException exc) {
            throw new ReviewedGATKException("Unable to close index file " + mFile, exc);
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
        final int bytesRequested = buffer.limit();
        if (bytesRequested == 0)
            return;

        try {

           //BufferedInputStream cannot read directly into a byte buffer, so we read into an array
            //and put the result into the bytebuffer after the if statement.

            // We have a rigid expectation here to read in exactly the number of bytes we've limited
            // our buffer to -- if there isn't enough data in the file, the index
            // must be truncated or otherwise corrupt:
            if(bytesRequested > fileLength - bufferedStream.position()){
                throw new UserException.MalformedFile(mFile, String.format("Premature end-of-file while reading BAM index file %s. " +
                        "It's likely that this file is truncated or corrupt -- " +
                        "Please try re-indexing the corresponding BAM file.",
                        mFile));
            }

            int bytesRead = bufferedStream.read(byteArray, 0, bytesRequested);

            // We have a rigid expectation here to read in exactly the number of bytes we've limited
            // our buffer to -- if we encounter EOF (-1), the index
            // must be truncated or otherwise corrupt:
            if (bytesRead <= 0) {
            throw new UserException.MalformedFile(mFile, String.format("Premature end-of-file while reading BAM index file %s. " +
                                                                       "It's likely that this file is truncated or corrupt -- " +
                                                                       "Please try re-indexing the corresponding BAM file.",
                                                                       mFile));
            }

            if(bytesRead != bytesRequested)
                throw new RuntimeException("Read amount different from requested amount. This should not happen.");

            buffer.put(byteArray, 0, bytesRequested);
        }
        catch(IOException ex) {
            throw new ReviewedGATKException("Index: unable to read bytes from index file " + mFile);
        }
    }


    /**
     * A reusable buffer for use by this index generator.
     * TODO: Should this be a SoftReference?
     */
    private ByteBuffer buffer = null;

    //BufferedStream don't read into ByteBuffers, so we need this temporary array
    private byte[] byteArray=null;
    private ByteBuffer getBuffer(final int size) {
        if(buffer == null || buffer.capacity() < size) {
            // Allocate a new byte buffer.  For now, make it indirect to make sure it winds up on the heap for easier debugging.
            buffer = ByteBuffer.allocate(size);
            byteArray = new byte[size];
            buffer.order(ByteOrder.LITTLE_ENDIAN);
        }
        buffer.clear();
        buffer.limit(size);
        return buffer;
    }

    private void skipBytes(final int count) {
        try {

            //try to skip forward the requested amount.
            long skipped =  bufferedStream.skip(count);

            if( skipped != count ) { //if not managed to skip the requested amount
                throw new ReviewedGATKException("Index: unable to reposition file channel of index file " + mFile);
            }
        }
        catch(IOException ex) {
            throw new ReviewedGATKException("Index: unable to reposition file channel of index file " + mFile);
        }
    }

    private void seek(final long position) {
        try {
            //to seek a new position, move the fileChannel, and reposition the bufferedStream
            bufferedStream.seek(position);
        }
        catch(IOException ex) {
            throw new ReviewedGATKException("Index: unable to reposition of file channel of index file " + mFile);
        }
    }

    /**
     * Retrieve the position from the current file channel.
     * @return position of the current file channel.
     */
    private long position() {
        try {
            return bufferedStream.position();
        }
        catch (IOException exc) {
            throw new ReviewedGATKException("Unable to read position from index file " + mFile, exc);
        }
    }
}
