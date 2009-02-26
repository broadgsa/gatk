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
 * Utility to count alignments (rather than gathering).
 */
public class CountAlignments {

    public static void main(String[] args)
        throws Exception {
        new CountAlignments().run(args);
    }

    private void usage() {
        System.out.println("Usage: CountAlignments ...");
        System.out.println("  -alignments <alignment-file> (- for stdin)");
        System.out.println("  -chromosome <chromosome>");
        System.out.println("  -start <start>");
        System.out.println("  -end <end>");
        System.out.println("  -bestAlignments");
        System.out.println("  -mismatchThreshold <n>");
        System.out.println("  -verbose");
        System.out.println("  -debug");
    }

    private boolean parseArguments(String[] args) {

        int argpos = 0;
        int argsleft = 0;

        while (argpos < args.length) {
            argsleft = args.length - argpos;
            String arg = args[argpos];
            if (arg.equals("-alignments") && argsleft > 1) {
                argpos++;
                mAlignmentFilePath = args[argpos++];
            } else if (arg.equals("-mismatchThreshold") && argsleft > 1) {
                argpos++;
                mMismatchThreshold = new Integer(args[argpos++]);
            } else if (arg.equals("-bestAlignments")) {
                argpos++;
                mReturnBestHits = true;
            } else if (arg.equals("-chromosome") && argsleft > 1) {
                argpos++;
                String chromosome = args[argpos++];
                mSequenceId = chromosomeToSequenceId(chromosome);
                if (mSequenceId < 0) {
                    System.out.println("Invalid chromosome: " + chromosome);
                    return false;
                }
            } else if (arg.equals("-start") && argsleft > 1) {
                argpos++;
                mStartPosition = new Integer(args[argpos++]);
            } else if (arg.equals("-end") && argsleft > 1) {
                argpos++;
                mEndPosition = new Integer(args[argpos++]);
            } else if (arg.equals("-verbose")) {
                argpos++;
                mVerbose = true;
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

        long[] counts = countAlignments(mAlignmentFilePath);
        String line = counts[0] + " " + counts[1];
        if (mAlignmentFilePath != null) {
            line = mAlignmentFilePath + " " + line;
        }
        System.out.println(line);
    }

    private long[] countAlignments(String path)
        throws IOException {
        long alignmentCount = 0;
        long baseCount = 0;
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
            if (mMismatchThreshold != null) {
                if (getAlignmentMismatches(alignment) > mMismatchThreshold) {
                    continue;
                }
            }
            if (mSequenceId != null) {
                if (alignment.getBSequenceId() != mSequenceId) {
                    continue;
                }
            }
            if (mStartPosition != null) {
                if (alignment.getBEnd() < mStartPosition) {
                    continue;
                }
            }
            if (mEndPosition != null) {
                if (alignment.getBStart() > mEndPosition) {
                    continue;
                }
            }
            alignmentCount++;
            baseCount += getBaseCount(alignment);
        }
        long[] result = { alignmentCount, baseCount };
        return result;
    }

    private Alignment getNextAlignment(LookAlignReader reader)
        throws IOException {
        if (!mReturnBestHits) {
            if (!reader.hasNext()) {
                return null;
            }
            return reader.next();
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
                return seed;
            }
            secondaryHits.add(seed);
            Alignment result = getUniqueBestAlignment(secondaryHits);
            if (result != null) {
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

    // Return the number of reference bases covered by this alignment.
    private int getBaseCount(Alignment alignment) {
        int count = 0;
        int[] blocks = alignment.getAlignmentBlocks();
        for (int i = 0; i < blocks.length; i += 3) {
            // int gap = blocks[i];
            int duration = blocks[i+1];
            // int mm = blocks[i+2];
            count += duration;
        }
        return count;
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

    private String mAlignmentFilePath = null;
    private boolean mReturnBestHits = false;
    private Integer mMismatchThreshold = null;
    private Integer mSequenceId = null;
    private Integer mStartPosition = null;
    private Integer mEndPosition = null;
    private Alignment mPendingAlignment = null;
}
