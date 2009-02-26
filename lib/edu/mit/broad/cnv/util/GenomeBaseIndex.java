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
 * Utility class for transforming between a linear base index
 * and a chromsome + position coordinate system.
 */
public class GenomeBaseIndex {

    private List<String> mSequenceNames = null;
    private int[] mLengths = null;
    private long[] mOffsets = null;

    private GenomeBaseIndex() {
    }

    public static GenomeBaseIndex read(File file)
        throws IOException {
        Reader reader = new BufferedReader(new FileReader(file));
        try {
            return read(reader);
        } finally {
            reader.close();
        }
    }

    // The input is just a list of space-delimited sequence name and length.
    public static GenomeBaseIndex read(Reader reader)
        throws IOException {
        List<String> sequenceNames = new ArrayList<String>();
        List<Integer> sequenceLengths = new ArrayList<Integer>();
        BufferedReader bufferedReader = new BufferedReader(reader);
        while (true) {
            String line = bufferedReader.readLine();
            if (line == null) {
                break;
            }
            String text = line.trim();
            if (text.length() == 0 || text.startsWith("#")) {
                continue;
            }
            String[] fields = text.split("\\s+");
            if (fields.length < 2) {
                throw new RuntimeException("Invalid input line: " + line);
            }
            int length = Integer.parseInt(fields[1]);
            if (length <= 0) {
                throw new RuntimeException("Invalid sequence length: " + length);
            }
            sequenceNames.add(fields[0]);
            sequenceLengths.add(length);
        }
        int count = sequenceLengths.size();
        int[] lengths = new int[count];
        long[] offsets = new long[count];
        long offset = 0;
        for (int i = 0; i < count; i++) {
            lengths[i] = sequenceLengths.get(i);
            offsets[i] = offset;
            offset += lengths[i];
        }
        GenomeBaseIndex result = new GenomeBaseIndex();
        result.mSequenceNames = sequenceNames;
        result.mLengths = lengths;
        result.mOffsets = offsets;
        return result;
    }

    public List<String> getSequenceNames() {
        return mSequenceNames;
    }

    public boolean contains(String seqName) {
        return (getSequenceIndex(seqName) >= 0);
    }

    public long getFirstIndex(String seqName) {
        int index = getSequenceIndex(seqName);
        if (index < 0) {
            return -1;
        }
        return mOffsets[index];
    }

    public long getLastIndex(String seqName) {
        int index = getSequenceIndex(seqName);
        if (index < 0) {
            return -1;
        }
        return (mOffsets[index] + mLengths[index] - 1);
    }

    public int getSequenceLength(String seqName) {
        int index = getSequenceIndex(seqName);
        if (index < 0) {
            return 0;
        }
        return mLengths[index];
    }

    public long getBaseIndex(String seqName, int position) {
        int index = getSequenceIndex(seqName);
        if (index < 0) {
            return -1;
        }
        if (position > mLengths[index]) {
            return -1;
        }
        if (position < 1) {
            // Zero or negative position means last base index
            position = mLengths[index];
        }
        return (mOffsets[index] + position - 1);
    }

    public String getSequenceName(long baseIndex) {
        int index = getSequenceIndex(baseIndex);
        if (index < 0) {
            return null;
        }
        return mSequenceNames.get(index);
    }

    public int getPosition(long baseIndex) {
        if (baseIndex < 0) {
            // Catch common sign-extension error when packing indexes as ints.
            throw new IllegalArgumentException("Invalid base index: " + baseIndex);
        }
        int index = getSequenceIndex(baseIndex);
        if (index < 0) {
            return 0;
        }
        long offset = mOffsets[index];
        long result = baseIndex - offset + 1;
        return (int) result;
    }

    // Same as getSequenceName, but treat the argument as an unsigned int.
    // This is useful for manipulating/storing indexes for the human
    // genome as 4-byte unsigned ints.
    public String getSequenceNameUnsigned(int baseIndex) {
        return getSequenceName(baseIndex & 0xFFFFFFFFL);
    }

    // Same as getPosition, but treat the argument as an unsigned int.
    // This is useful for manipulating/storing indexes for the human
    // genome as 4-byte unsigned ints.
    public int getPositionUnsigned(int baseIndex) {
        return getPosition(baseIndex & 0xFFFFFFFFL);
    }

    private int getSequenceIndex(String seqName) {
        return mSequenceNames.indexOf(seqName);
    }

    private int getSequenceIndex(long baseIndex) {
        long offset = 0;
        if (baseIndex < 0) {
            return -1;
        }
        for (int i = 0; i < mLengths.length; i++) {
            int length = mLengths[i];
            if (offset + length > baseIndex) {
                return i;
            }
            offset += length;
        }
        return -1;
    }
}
