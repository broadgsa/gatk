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
package edu.mit.broad.cnv;

import edu.mit.broad.arachne.Alignment;
import edu.mit.broad.arachne.LookAlignReader;

import java.io.*;
import java.util.*;


/**
 * Utility class to do data reduction on CNV data.
 */
public class AnalyzeCnvs {

    public static void main(String[] args)
        throws Exception {
        new AnalyzeCnvs().run(args);
    }

    private void usage() {
        System.out.println("Usage: AnalyzeCnvs ...");
        System.out.println("  -action <action>");
        System.out.println("  -alignments <alignment-file> or -");
        System.out.println("  -alignmentList <alignment-fofn>");
        System.out.println("  -chromosome <chrN>");
        System.out.println("  -start <start-coordinate>");
        System.out.println("  -end <end-coordinate>");
        System.out.println("  -bestAlignments");
        System.out.println("  -mismatchThreshold <n>");
        System.out.println("  -binsize <n>");
        System.out.println("  -output <coverage|all>");
        System.out.println("  -verbose");
        System.out.println("  -debug");
    }

    private boolean parseArguments(String[] args) {

        int argpos = 0;
        int argsleft = 0;

        while (argpos < args.length) {
            argsleft = args.length - argpos;
            String arg = args[argpos];
            if (arg.equals("-action") && argsleft > 1) {
                argpos++;
                mAction = args[argpos++];
            } else if (arg.equals("-alignments") && argsleft > 1) {
                argpos++;
                mAlignmentFilePath = args[argpos++];
            } else if (arg.equals("-alignmentList") && argsleft > 1) {
                argpos++;
                mAlignmentListFilePath = args[argpos++];
            } else if (arg.equals("-chromosome") && argsleft > 1) {
                argpos++;
                mChromosome = args[argpos++];
            } else if (arg.equals("-start") && argsleft > 1) {
                argpos++;
                mStartPosition = new Integer(args[argpos++]);
            } else if (arg.equals("-end") && argsleft > 1) {
                argpos++;
                mEndPosition = new Integer(args[argpos++]);
            } else if (arg.equals("-verbose")) {
                argpos++;
                mVerbose = true;
            } else if (arg.equals("-mismatchThreshold") && argsleft > 1) {
                argpos++;
                mMismatchThreshold = new Integer(args[argpos++]);
            } else if (arg.equals("-bestAlignments")) {
                argpos++;
                mReturnBestHits = true;
            } else if (arg.equals("-binsize") && argsleft > 1) {
                argpos++;
                mBinSize = Integer.parseInt(args[argpos++]);
            } else if (arg.equals("-output") && argsleft > 1) {
                argpos++;
                mOutputColumns = args[argpos++];
            } else if (arg.equals("-debug")) {
                argpos++;
                mDebug = true;
            } else if (arg.startsWith("-")) {
                usage();
                return false;
            } else {
                break;
            }
        }

        argsleft = args.length - argpos;
        if (argsleft != 0) {
            usage();
            return false;
        }

        return true;
    }

    private void run(String[] args)
        throws Exception {

        if (!parseArguments(args)) {
            System.exit(1);
        }

        if (mAction == null) {
            mAction = "alignmentCoverage";
        }

        if (mAction.equals("alignmentCoverage")) {
            mainAlignmentCoverage();
        } else {
            System.out.println("Unknown action: " + mAction);
            usage();
            System.exit(1);
        }
    }

    private void mainAlignmentCoverage()
        throws IOException {

        if (mStartPosition == null || mEndPosition == null) {
            usage();
            System.exit(1);
        } else if (mStartPosition <= 0 || mEndPosition <= 0 || mStartPosition > mEndPosition) {
            System.out.println("Invalid start/end positions: " + mStartPosition + " " + mEndPosition);
            usage();
            System.exit(1);
        }

        mSequenceId = chromosomeToSequenceId(mChromosome);
        if (mSequenceId < 0) {
            System.out.println("Invalid chromosome: " + mChromosome);
            usage();
            System.exit(1);
        }

        if (mBinSize <= 0) {
            System.out.println("Invalid bin size: " + mBinSize);
            usage();
            System.exit(1);
        }

        runAlignmentCoverage();
    }

    private void runAlignmentCoverage()
        throws IOException {

        int length = (mEndPosition - mStartPosition + 1);
        if (length <= 0) {
            throw new RuntimeException("Invalid start/end positions");
        }

        int binSize = mBinSize;
        int binCount = (length + binSize - 1) / binSize;
        int[] readStarts = new int[binCount];
        int[] readDepths = new int[binCount];
        List<String> alignmentFiles = getAlignmentFiles();
        for (String path : alignmentFiles) {
            processAlignmentFile(path, readStarts, readDepths);
        }
        printStats(readStarts, readDepths);
    }

    private List<String> getAlignmentFiles()
        throws IOException {
        List<String> fileList = new ArrayList<String>();
        if (mAlignmentListFilePath != null) {
            LineNumberReader reader = new LineNumberReader(new FileReader(mAlignmentListFilePath));
            while (true) {
                String line = reader.readLine();
                if (line == null) {
                    reader.close();
                    break;
                }
                String path = line.trim();
                if (path.length() == 0 || path.startsWith("#")) {
                    continue;
                }
                fileList.add(path);
            }
        } else if (mAlignmentFilePath != null) {
            fileList.add(mAlignmentFilePath);
        }
        return fileList;
    }

    private void processAlignmentFile(String path, int[] readStarts, int[] readDepths)
        throws IOException {

        LookAlignReader reader = null;
        if (path == null || path.equals("-")) {
            reader = new LookAlignReader(new InputStreamReader(System.in));
        } else {
            reader = new LookAlignReader(new File(path));
        }

        while (true) {
            Alignment alignment = getNextAlignment(reader);
            if (alignment == null) {
                reader.close();
                break;
            }
            processAlignment(alignment, readStarts, readDepths);
        }
    }

    private void processAlignment(Alignment alignment,
                                  int[] readStarts,
                                  int[] readDepths) {

        if (readStarts != null) {
            int baseOffset = alignment.getBStart() - mStartPosition;
            int binIndex = baseOffset / mBinSize;
            if (binIndex >= 0 && binIndex < readStarts.length) {
                readStarts[binIndex]++;
            }
        }

        if (readDepths != null) {
            int baseOffset = alignment.getBStart() - mStartPosition;
            int[] alignmentBlocks = alignment.getAlignmentBlocks();
            for (int i = 0; i < alignmentBlocks.length; i += 3) {
                int gap = alignmentBlocks[i];
                int duration = alignmentBlocks[i+1];
                if (gap > 0) {
                    // Gap in B sequence (genome)
                    // Negative gaps are gaps in A sequence (read)
                    baseOffset += gap;
                }
                for (int j = 0; j < duration; j++) {
                    int binIndex = baseOffset / mBinSize;
                    if (binIndex >= 0 && binIndex < readDepths.length) {
                        readDepths[binIndex]++;
                    }
                    baseOffset++;
                }
            }
        }
    }

    private Alignment getNextAlignment(LookAlignReader reader)
        throws IOException {

        if (!mReturnBestHits) {
            while (reader.hasNext()) {
                Alignment alignment = reader.next();
                if (passesAlignmentFilters(alignment)) {
                    return alignment;
                }
            }
            return null;
        }

        while (true) {
            Alignment seed = mPendingAlignment;
            mPendingAlignment = null;
            if (seed == null && reader.hasNext()) {
                seed = reader.next();
            }
            if (seed == null) {
                return null;
            }
            List<Alignment> secondaryHits = null;
            while (reader.hasNext()) {
                Alignment alignment = reader.next();
                if (alignment.getASequenceId() != seed.getASequenceId()) {
                    if (alignment.getASequenceId() < seed.getASequenceId()) {
                        throw new RuntimeException("Alignments not sorted by A sequence: " + alignment.format());
                    }
                    mPendingAlignment = alignment;
                    break;
                }
                if (secondaryHits == null) {
                    secondaryHits = new ArrayList<Alignment>();
                }
                secondaryHits.add(alignment);
            }
            if (secondaryHits == null) {
                if (!passesAlignmentFilters(seed)) {
                    continue;
                }
                return seed;
            }
            secondaryHits.add(seed);
            Alignment result = getUniqueBestAlignment(secondaryHits);
            if (result != null && passesAlignmentFilters(result)) {
                return result;
            }
        }
    }

    private Alignment getUniqueBestAlignment(List<Alignment> alignments) {
        int bestMismatches = 0;
        List<Alignment> best = new ArrayList<Alignment>();
        for (Alignment a : alignments) {
            int mismatches = getAlignmentMismatches(a);
            if (best.isEmpty()) {
                best.add(a);
                bestMismatches = mismatches;
            }
            if (mismatches == bestMismatches) {
                best.add(a);
            } else if (mismatches < bestMismatches) {
                best.clear();
                best.add(a);
                bestMismatches = mismatches;
            }
        }
        if (best.size() != 1) {
            return null;
        }
        return best.get(0);
    }

    private boolean passesAlignmentFilters(Alignment alignment) {

        if (mMismatchThreshold != null) {
            if (getAlignmentMismatches(alignment) > mMismatchThreshold) {
                return false;
            }
        }

        if (mSequenceId != null) {
            if (alignment.getBSequenceId() != mSequenceId) {
                return false;
            }
        }

        if (mStartPosition != null) {
            if (alignment.getBEnd() < mStartPosition) {
                return false;
            }
        }

        if (mEndPosition != null) {
            if (alignment.getBStart() > mEndPosition) {
                return false;
            }
        }

        return true;
    }

    private int getAlignmentMismatches(Alignment alignment) {
        int mismatches = 0;
        int[] blocks = alignment.getAlignmentBlocks();
        for (int i = 0; i < blocks.length; i += 3) {
            int gap = blocks[i];
            int duration = blocks[i+1];
            int mm = blocks[i+2];
            if (mm > duration) {
                throw new RuntimeException("Invalid alignment? : " + alignment.format());
            }
            mismatches += Math.abs(gap);
            mismatches += mm;
        }
        return mismatches;
    }

    private void printStats(int[] readStarts, int[] readDepths) {
        if (mOutputColumns != null && mOutputColumns.equals("coverage")) {
            // No headers, just coverage
            for (int i = 0; i < readDepths.length; i++) {
                String line = "";
                if (mBinSize == 1) {
                    line += readDepths[i];
                } else {
                    line += (readDepths[i] / (double) mBinSize);
                }
                System.out.println(line);
            }
        } else {
            System.out.println("Position" + "\t" + "Starts" + "\t" + "Coverage");
            for (int i = 0; i < readDepths.length; i++) {
                String line = "";
                int position = mStartPosition + i*mBinSize;
                line += position + "\t" + readStarts[i] + "\t";
                if (mBinSize == 1) {
                    line += readDepths[i];
                } else {
                    line += (readDepths[i] / (double) mBinSize);
                }
                System.out.println(line);
            }
        }
    }

    private int chromosomeToSequenceId(String text) {
        if (text == null || text.length() == 0) {
            return -1;
        }
        if (text.matches("\\d+")) {
            return Integer.parseInt(text);
        }
        if (text.startsWith("chr") && text.length() > 3) {
            text = text.substring(3);
        }
        if (text.matches("\\d+") && !text.startsWith("0")) {
            return Integer.parseInt(text);
        }
        if (text.equals("M")) {
            return 0;
        } else if (text.equals("X")) {
            return 23;
        } else if (text.equals("Y")) {
            return 24;
        } else {
            return -1;
        }
    }

    private boolean mDebug = false;
    private boolean mVerbose = false;

    private String mAction = null;
    private String mAlignmentFilePath = null;
    private String mAlignmentListFilePath = null;
    private String mChromosome = null;
    private Integer mStartPosition = null;
    private Integer mEndPosition = null;
    private Integer mSequenceId = null;
    private boolean mReturnBestHits = false;
    private Integer mMismatchThreshold = null;
    private int mBinSize = 1;
    private String mOutputColumns = null;
    private Alignment mPendingAlignment = null;
}
