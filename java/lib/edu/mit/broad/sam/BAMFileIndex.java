/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2008 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package edu.mit.broad.sam;


import edu.mit.broad.sam.util.RuntimeEOFException;
import edu.mit.broad.sam.util.RuntimeIOException;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.util.*;

/**
 * Internal class for reading BAM file indexes.
 */
class BAMFileIndex
{
    private static final int MAX_BINS = 37450; // =(8^6-1)/7+1
    private static final int BAM_LIDX_SHIFT = 16;

    private File mFile = null;
    private FileInputStream mFileStream = null;
    private MappedByteBuffer mFileBuffer = null;


    BAMFileIndex(final File file) {
        mFile = file;
    }

    void close() {
        closeFileStream();
    }

    long[] getSearchBins(int referenceIndex, int startPos, int endPos) {

        openIndex();
        seek(4);

        int sequenceCount = readInteger();
        // System.out.println("# Sequence count: " + sequenceCount);
        if (referenceIndex >= sequenceCount) {
            return null;
        }

        BitSet regionBins = regionToBins(startPos, endPos);
        if (regionBins == null) {
            return null;
        }

        for (int i = 0; i < referenceIndex; i++) {
            // System.out.println("# Sequence TID: " + i);
            int nBins = readInteger();
            // System.out.println("# nBins: " + nBins);
            for (int j = 0; j < nBins; j++) {
                int bin = readInteger();
                int nChunks = readInteger();
                // System.out.println("# bin[" + j + "] = " + bin + ", nChunks = " + nChunks);
                skipBytes(16 * nChunks);
            }
            int nLinearBins = readInteger();
            // System.out.println("# nLinearBins: " + nLinearBins);
            skipBytes(8 * nLinearBins);
        }

        // System.out.println("# Sequence target TID: " + referenceIndex);
        int nIndexBins = readInteger();
        // System.out.println("# nBins: " + nIndexBins);
        if (nIndexBins == 0) {
            return null;
        }

        List<Chunk> chunkList = new ArrayList<Chunk>();
        for (int i = 0; i < nIndexBins; i++) {
            int indexBin = readInteger();
            int nChunks = readInteger();
            // System.out.println("# bin[" + i + "] = " + indexBin + ", nChunks = " + nChunks);
            if (regionBins.get(indexBin)) {
                for (int ci = 0; ci < nChunks; ci++) {
                    long chunkBegin = readLong();
                    long chunkEnd = readLong();
                    chunkList.add(new Chunk(chunkBegin, chunkEnd));
                }
            } else {
                skipBytes(16 * nChunks);
            }
        }

        if (chunkList.isEmpty()) {
            return null;
        }

        int start = (startPos <= 0) ? 0 : startPos-1;
        int regionLinearBin = start >> BAM_LIDX_SHIFT;
        int nLinearBins = readInteger();
        // System.out.println("# nLinearBins: " + nLinearBins);
        // System.out.println("# regionLinearBin: " + regionLinearBin);
        long minimumOffset = 0;
        if (regionLinearBin < nLinearBins) {
            skipBytes(8 * regionLinearBin);
            minimumOffset = readLong();
        }
        chunkList = optimizeChunkList(chunkList, minimumOffset);
        return convertToArray(chunkList);
    }

    private List<Chunk> optimizeChunkList(List<Chunk> chunkList, long minimumOffset) {
        Chunk lastChunk = null;
        Collections.sort(chunkList);
        List<Chunk> result = new ArrayList<Chunk>();
        for (Chunk chunk : chunkList) {
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
            long lastFileBlock = getFileBlock(lastChunk.getChunkEnd());
            long chunkFileBlock = getFileBlock(chunk.getChunkStart());
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

    private long[] convertToArray(List<Chunk> chunkList) {
        int count = chunkList.size() * 2;
        if (count == 0) {
            return null;
        }
        int index = 0;
        long[] result = new long[count];
        for (Chunk chunk : chunkList) {
            result[index++] = chunk.getChunkStart();
            result[index++] = chunk.getChunkEnd();
        }
        return result;
    }

    private BitSet regionToBins(int startPos, int endPos) {
        int maxPos = 0x1FFFFFFF;
        int start = (startPos <= 0) ? 0 : (startPos-1) & maxPos;
        int end = (endPos <= 0) ? maxPos : (endPos-1) & maxPos;
        if (start > end) {
            return null;
        }
        int k;
        BitSet bitSet = new BitSet(MAX_BINS);
        bitSet.set(0);
        for (k =    1 + (start>>26); k <=    1 + (end>>26); ++k) bitSet.set(k);
        for (k =    9 + (start>>23); k <=    9 + (end>>23); ++k) bitSet.set(k);
        for (k =   73 + (start>>20); k <=   73 + (end>>20); ++k) bitSet.set(k);
        for (k =  585 + (start>>17); k <=  585 + (end>>17); ++k) bitSet.set(k);
        for (k = 4681 + (start>>14); k <= 4681 + (end>>14); ++k) bitSet.set(k);
        return bitSet;
    }

    private long getFileBlock(long bgzfOffset) {
        return ((bgzfOffset >> 16L) & 0xFFFFFFFFFFFFL);
    }

    private void openIndex() {
        if (mFileBuffer != null) {
            return;
        }
        openFileStream();
        seek(0);
        byte[] buffer = new byte[4];
        readBytes(buffer);
        if (!Arrays.equals(buffer, BAMFileConstants.BAM_INDEX_MAGIC)) {
            closeFileStream();
            throw new RuntimeException("Invalid file header in BAM index " + mFile +
                                       ": " + new String(buffer));
        }
    }

    private void readBytes(byte[] buffer) {
        mFileBuffer.get(buffer);
    }

    private int readInteger() {
        return mFileBuffer.getInt();
    }

    private long readLong() {
        return mFileBuffer.getLong();
    }

    private void skipBytes(int count) {
        mFileBuffer.position(mFileBuffer.position() + count);
    }

    private void seek(int position) {
        mFileBuffer.position(position);
    }

    private void openFileStream() {
        if (mFileStream != null) {
            return;
        }
        try {
            mFileStream = new FileInputStream(mFile);
            FileChannel channel = mFileStream.getChannel();
            mFileBuffer = channel.map(FileChannel.MapMode.READ_ONLY, 0L, channel.size());
            mFileBuffer.order(ByteOrder.LITTLE_ENDIAN);
        } catch (IOException exc) {
            throw new RuntimeIOException(exc.getMessage(), exc);
        }
    }

    private void closeFileStream() {
        if (mFileStream == null) {
            return;
        }
        try {
            mFileStream.close();
        } catch (IOException exc) {
            throw new RuntimeIOException(exc.getMessage(), exc);
        }
        mFileStream = null;
        mFileBuffer = null;
    }

    private static class Chunk
        implements Comparable<Chunk> {

        private long mChunkStart;
        private long mChunkEnd;

        Chunk(long start, long end) {
            mChunkStart = start;
            mChunkEnd = end;
        }

        long getChunkStart() {
            return mChunkStart;
        }

        void setChunkStart(long value) {
            mChunkStart = value;
        }

        long getChunkEnd() {
            return mChunkEnd;
        }

        void setChunkEnd(long value) {
            mChunkEnd = value;
        }

        public int compareTo(Chunk chunk) {
            int result = Long.signum(mChunkStart - chunk.mChunkStart);
            if (result == 0) {
                result = Long.signum(mChunkEnd - chunk.mChunkEnd);
            }
            return result;
        }
    }
}
