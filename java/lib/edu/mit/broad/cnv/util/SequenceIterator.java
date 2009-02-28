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
 * Utility class for iterating over fasta files.
 * Also maintains an unsigned base index over the file set.
 */
public class SequenceIterator
{
    private List<File> mInputFiles = null;
    private int mInputFileIndex = 0;
    private int mBaseIndex = -1;
    private LineNumberReader mCurrentReader = null;
    private String mNextSequence = null;
    private String mLineBuffer = null;
    private int mLineBufferIndex = 0;

    public SequenceIterator(File inputFile) {
        mInputFiles = new ArrayList<File>();
        mInputFiles.add(inputFile);
    }

    public SequenceIterator(List<File> inputFiles) {
        mInputFiles = inputFiles;
    }

    public void close() {
        if (mCurrentReader != null) {
            try {
                mCurrentReader.close();
            } catch (IOException exc) {
                throw new RuntimeException("Error closing reader: " + exc.getMessage(),
                                           exc);
            }
        }
        mCurrentReader = null;
        mInputFiles = null;
        mInputFileIndex = 0;
        mBaseIndex = -1;
        mNextSequence = null;
        mLineBuffer = null;
        mLineBufferIndex = 0;
    }

    public String getNextSequence()
        throws IOException {

        while (mNextSequence == null) {
            if (mLineBuffer != null) {
                incrementBaseIndex(mLineBuffer.length() - mLineBufferIndex);
                mLineBuffer = null;
                mLineBufferIndex = 0;
            }
            if (mCurrentReader == null) {
                mCurrentReader = getNextReader();
                if (mCurrentReader == null) {
                    return null;
                }
            }
            String line = mCurrentReader.readLine();
            if (line == null) {
                mCurrentReader.close();
                mCurrentReader = null;
                continue;
            }
            if (line.startsWith(">")) {
                String[] tokens = line.substring(1).trim().split("\\s+");
                mNextSequence = tokens[0];
            } else {
                incrementBaseIndex(line.length());
            }
        }
        String result = mNextSequence;
        mNextSequence = null;
        return result;
    }

    public char getNextBase()
        throws IOException {

        if (mLineBuffer == null || mLineBufferIndex >= mLineBuffer.length()) {
            if (mCurrentReader == null) {
                return 0;
            }
            if (mNextSequence != null) {
                return 0;
            }
            String line = mCurrentReader.readLine();
            if (line == null) {
                mLineBuffer = null;
                mLineBufferIndex = 0;
                mCurrentReader.close();
                mCurrentReader = null;
                return 0;
            }
            if (line.startsWith(">")) {
                String[] tokens = line.substring(1).trim().split("\\s+");
                mNextSequence = tokens[0];
                mLineBuffer = null;
                mLineBufferIndex = 0;
                return 0;
            }
            mLineBuffer = line.toUpperCase();
            mLineBufferIndex = 0;
        }
        char result = mLineBuffer.charAt(mLineBufferIndex++);
        incrementBaseIndex(1);
        return result;
    }

    public int getBaseIndex() {
        return mBaseIndex;
    }

    private LineNumberReader getNextReader()
        throws IOException {
        if (mInputFileIndex >= mInputFiles.size()) {
            return null;
        }
        File file = mInputFiles.get(mInputFileIndex++);
        return new LineNumberReader(new FileReader(file));
    }

    private void incrementBaseIndex(int amount) {
        if (mBaseIndex < -1 && (mBaseIndex + amount) >= -1) {
            throw new RuntimeException("Base index: 32-bit overflow");
        }
        mBaseIndex += amount;
    }
}

