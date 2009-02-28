/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package edu.mit.broad.cnv.util;


import java.io.*;
import java.util.*;


/**
 * Utility class for transforming between a chromsome + position
 * coordinate system and a binned coordinate system where each
 * chromosome (separately) is divided into fixed sized bins,
 * ragged on the right/upper end.
 */
public class GenomeBinIndex {

    private int mBinSize;
    private List<String> mSequenceNames;
    private int[] mSequenceLengths;
    private int[] mBinOffsets;

    public GenomeBinIndex(GenomeBaseIndex gbi, int binSize) {
        if (binSize <= 0) {
            throw new IllegalArgumentException("Illegal bin size: " + binSize);
        }
        mBinSize = binSize;
        mSequenceNames = new ArrayList<String>(gbi.getSequenceNames());
        int count = mSequenceNames.size();
        mSequenceLengths = new int[count];
        mBinOffsets = new int[count];
        long binOffset = 0;  // long to detect overflow
        for (int i = 0; i < count; i++) {
            int length = gbi.getSequenceLength(mSequenceNames.get(i));
            int binCount = (length + binSize - 1) / binSize;
            mSequenceLengths[i] = length;
            mBinOffsets[i] = (int) binOffset;
            binOffset += binCount;
        }
        if (binOffset > Integer.MAX_VALUE) {
            // Check for integer overflow.
            // This will happen, e.g., with the human genome and a bin size of 1.
            throw new RuntimeException("Binsize too small: " + binSize);
        }
    }

    public int getBinSize() {
        return mBinSize;
    }

    public int getBinIndex(String seqName, int position) {
        int index = getSequenceIndex(seqName);
        if (index < 0) {
            return -1;
        }
        if (position > mSequenceLengths[index]) {
            return -1;
        }
        if (position < 1) {
            position = mSequenceLengths[index];
        }
        int bin = (position - 1) / mBinSize;
        return (mBinOffsets[index] + bin);
    }

    public String getSequenceName(int binIndex) {
        int index = getSequenceIndex(binIndex);
        if (index < 0) {
            return null;
        }
        return mSequenceNames.get(index);
    }

    public int getStartPosition(int binIndex) {
        int index = getSequenceIndex(binIndex);
        if (index < 0) {
            return -1;
        }
        int bin = binIndex - mBinOffsets[index];
        return (bin * mBinSize + 1);
    }

    public int getEndPosition(int binIndex) {
        int index = getSequenceIndex(binIndex);
        if (index < 0) {
            return -1;
        }
        int bin = binIndex - mBinOffsets[index];
        int position = (bin+1) * mBinSize;
        position = Math.min(position, mSequenceLengths[index]);
        return position;
    }

    public List<String> getSequenceNames() {
        return mSequenceNames;
    }

    public int getFirstBin(String seqName) {
        return getBinIndex(seqName, 1);
    }

    public int getLastBin(String seqName) {
        return getBinIndex(seqName, 0);
    }

    public int getBinCount() {
        if (mBinOffsets.length == 0) {
            return 0;
        }
        int lastIndex = mBinOffsets.length - 1;
        int count = mBinOffsets[lastIndex];
        count += (mSequenceLengths[lastIndex] + mBinSize - 1) / mBinSize;
        return count;
    }

    public int getBinCount(String seqName) {
        int index = getSequenceIndex(seqName);
        if (index < 0) {
            return -1;
        }
        return ((mSequenceLengths[index] + mBinSize - 1) / mBinSize);
    }

    public int getSequenceLength(String seqName) {
        int index = getSequenceIndex(seqName);
        if (index < 0) {
            return 0;
        }
        return mSequenceLengths[index];
    }

    private int getSequenceIndex(String seqName) {
        for (int i = 0; i < mSequenceNames.size(); i++) {
            if (mSequenceNames.get(i).equals(seqName)) {
                return i;
            }
        }
        return -1;
    }

    private int getSequenceIndex(int binIndex) {
        if (binIndex < 0) {
            return -1;
        }
        for (int i = 1; i < mBinOffsets.length; i++) {
            if (mBinOffsets[i] > binIndex) {
                return i-1;
            }
        }
        int lastIndex = mBinOffsets.length-1;
        int lastBinIndex = mBinOffsets[lastIndex];
        lastBinIndex += (mSequenceLengths[lastIndex] + mBinSize - 1) / mBinSize;
        if (binIndex <= lastBinIndex) {
            return lastIndex;
        }
        return -1;
    }
}

