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
package edu.mit.broad.cnv.kmer;


import edu.mit.broad.dcp.DistributedAlgorithm;
import edu.mit.broad.cnv.util.SequenceIterator;

import java.io.*;
import java.util.*;


/**
 * Distributed algorithm for counting unique kmers.
 */
public class DistributedKMerCounter
    extends DistributedAlgorithm
{
    private boolean mDebug = false;
    private boolean mVerbose = false;
    private int mK = 0;
    private List<File> mInputFiles = null;
    private List<String> mSequenceList = null;
    private List<Integer> mSequenceOffsetList = null;


    public DistributedKMerCounter() {
    }

    public boolean getDebug() {
        return mDebug;
    }

    public void setDebug(boolean value) {
        mDebug = value;
    }

    public boolean getVerbose() {
        return mVerbose;
    }

    public void setVerbose(boolean value) {
        mVerbose = value;
    }

    public int getK() {
        return mK;
    }

    public void setK(int value) {
        mK = value;
    }

    public List<File> getInputFiles() {
        return mInputFiles;
    }

    public void setInputFiles(List<File> value) {
        mInputFiles = value;
    }

    public void run()
        throws Exception {
        super.run();
        finish();
    }

    protected void init()
        throws Exception {
        if (getWorkerId() == MASTER) {
            initMaster();
        } else {
            initWorker();
        }
    }

    private void initMaster()
        throws IOException {
        // Tasks to be amortized
        report("Scanning sequences ...");
        scanSequences();
        report("Scan complete.");
    }

    private void initWorker() {
        // Tasks to be amortized
    }

    protected void start() {
        // scan genome, divide into chromosomes and optionally segments, distribute calls
    }

    private void finish() {
        // merge individual files, write out final results
    }

    private void scanSequences()
        throws IOException {
        List<String> sequenceList = new ArrayList<String>();
        List<Integer> sequenceOffsetList = new ArrayList<Integer>();
        SequenceIterator seqIterator = new SequenceIterator(getInputFiles());
        while (true) {
            String seqName = seqIterator.getNextSequence();
            if (seqName == null) {
                break;
            }
            int baseIndex = seqIterator.getBaseIndex() + 1;
            sequenceList.add(seqName);
            sequenceOffsetList.add(baseIndex);
        }
        mSequenceList = sequenceList;
        mSequenceOffsetList = sequenceOffsetList;
    }

    // Currently not used
    private void loadGenomeOffsets(File file)
        throws IOException {
        List<String> sequenceList = new ArrayList<String>();
        List<Integer> sequenceOffsetList = new ArrayList<Integer>();
        int baseIndex = 0;
        LineNumberReader reader = new LineNumberReader(new FileReader(file));
        while (true) {
            String line = reader.readLine();
            if (line == null) {
                break;
            }
            String text = line.trim();
            if (text.length() == 0 || text.startsWith("#")) {
                continue;
            }
            String[] fields = text.split("\\s+");
            if (fields.length != 2) {
                throw new RuntimeException("Invalid input line: " + line);
            }
            int length = Integer.parseInt(fields[1]);
            sequenceList.add(fields[0]);
            sequenceOffsetList.add(baseIndex);
            baseIndex += length;
        }
        mSequenceList = sequenceList;
        mSequenceOffsetList = sequenceOffsetList;
    }
}
