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
 * Utility program to gather CNV alignments from LookAlign files in an I/O efficient manner.
 */
public class GatherAlignments {

    public static void main(String[] args)
        throws Exception {
        new GatherAlignments().run(args);
    }

    private void usage() {
        System.out.println("Usage: GatherAlignments ...");
        System.out.println("  -cnpList <cnp-file>");
        System.out.println("  -sampleId <sample-id>");
        System.out.println("  -inputFileList <fofn>");
        System.out.println("  -outputDirectory <dir>");
        System.out.println("  -padding <n-bases>");
        System.out.println("  -bestAlignments");
        System.out.println("  -verbose");
        System.out.println("  -debug");
    }

    private boolean parseArguments(String[] args) {

        int argpos = 0;
        int argsleft = 0;

        while (argpos < args.length) {
            argsleft = args.length - argpos;
            String arg = args[argpos];
            if (arg.equals("-cnpList") && argsleft > 1) {
                argpos++;
                mCnpListPath = args[argpos++];
            } else if (arg.equals("-sampleId") && argsleft > 1) {
                argpos++;
                mSampleId = args[argpos++];
            } else if (arg.equals("-inputFileList") && argsleft > 1) {
                argpos++;
                mInputFileListPath = args[argpos++];
            } else if (arg.equals("-outputDirectory") && argsleft > 1) {
                argpos++;
                mOutputDirectory = args[argpos++];
            } else if (arg.equals("-padding") && argsleft > 1) {
                argpos++;
                mCnpRegionPadding = Integer.parseInt(args[argpos++]);
            } else if (arg.equals("-bestAlignments")) {
                argpos++;
                mReturnBestHits = true;
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

        List<File> mInputFileList = parseInputFiles(mInputFileListPath);
        Map<Integer, List<CnpRegion>> mCnpMap = parseCnpFile(mCnpListPath);
        for (File inputFile : mInputFileList) {
            scanInputFile(inputFile, mCnpMap);
        }
    }

    private List<File> parseInputFiles(String path)
        throws IOException {
        List<File> fileList = new ArrayList<File>();
        LineNumberReader reader = new LineNumberReader(new FileReader(path));
        while (true) {
            String line = reader.readLine();
            if (line == null) {
                reader.close();
                break;
            }
            line = line.trim();
            if (line.length() == 0 || line.startsWith("#")) {
                continue;
            }
            String[] fields = line.split("\\s+");
            fileList.add(new File(fields[0]));
        }
        return fileList;
    }

    private Map<Integer, List<CnpRegion>> parseCnpFile(String path)
        throws IOException {
        Map<Integer, List<CnpRegion>> cnpMap = new HashMap<Integer, List<CnpRegion>>();
        LineNumberReader reader = new LineNumberReader(new FileReader(path));
        while (true) {
            String line = reader.readLine();
            if (line == null) {
                reader.close();
                break;
            }
            line = line.trim();
            if (line.length() == 0 || line.startsWith("#")) {
                continue;
            }
            String[] fields = line.split("\\s+");
            if (fields.length != 4) {
                throw new RuntimeException("Invalid CNP line: " + line);
            }
            if (fields[0].equalsIgnoreCase("CNPID")) {
                continue;
            }
            String cnpId = fields[0];
            String chromosome = fields[1];
            int start = Integer.parseInt(fields[2].replaceAll(",", ""));
            int end = Integer.parseInt(fields[3].replaceAll(",", ""));
            int sequenceId = chromosomeToSequenceId(chromosome);
            if (sequenceId < 0) {
                throw new RuntimeException("Unrecognized chromosome: " + chromosome);
            }
            if (mCnpRegionPadding > 0) {
                start = Math.max(1, start - mCnpRegionPadding);
                end = end + mCnpRegionPadding;
            }
            CnpRegion cnp = new CnpRegion(cnpId, sequenceId, start, end);
            List<CnpRegion> cnpList = cnpMap.get(sequenceId);
            if (cnpList == null) {
                cnpList = new ArrayList<CnpRegion>();
                cnpMap.put(sequenceId, cnpList);
            }
            cnpList.add(cnp);
        }
        return cnpMap;
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

    private void scanInputFile(File inputFile,
                               Map<Integer, List<CnpRegion>> cnpMap)
        throws IOException {
        LookAlignReader reader = new LookAlignReader(inputFile);
        while (true) {
            Alignment alignment = getNextAlignment(reader);
            if (alignment == null) {
                reader.close();
                break;
            }
            List<CnpRegion> cnpList = cnpMap.get(alignment.getBSequenceId());
            if (cnpList == null) {
                continue;
            }
            for (CnpRegion cnp : cnpList) {
                if (overlaps(cnp, alignment)) {
                    saveCnpAlignment(cnp, alignment, inputFile);
                }
            }
        }
        flushCnpAlignments(inputFile);
    }

    private Alignment getNextAlignment(LookAlignReader reader)
        throws IOException {
        if (!mReturnBestHits) {
            if (reader.hasNext()) {
                return reader.next();
            } else {
                return null;
            }
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

    private boolean overlaps(CnpRegion cnp, Alignment alignment) {
        return (cnp.getSequenceId() == alignment.getBSequenceId() &&
                cnp.getStart() <= alignment.getBEnd() &&
                cnp.getEnd() >= alignment.getBStart());
    }

    private void saveCnpAlignment(CnpRegion cnp, Alignment alignment, File inputFile)
        throws IOException {
        if (mCnpAlignmentCount > mCnpAlignmentLimit) {
            flushCnpAlignments(inputFile);
        }
        String cnpId = cnp.getCnpId();
        List<Alignment> alignmentList = mCnpAlignmentMap.get(cnpId);
        if (alignmentList == null) {
            alignmentList = new ArrayList<Alignment>();
            mCnpAlignmentMap.put(cnpId, alignmentList);
        }
        alignmentList.add(alignment);
        mCnpAlignmentCount++;
    }

    private void flushCnpAlignments(File inputFile)
        throws IOException {
        while (!mCnpAlignmentMap.isEmpty()) {
            String cnpId = mCnpAlignmentMap.keySet().iterator().next();
            List<Alignment> alignmentList = mCnpAlignmentMap.get(cnpId);
            writeAlignments(cnpId, mSampleId, alignmentList, inputFile);
            mCnpAlignmentMap.remove(cnpId);
            mCnpAlignmentCount -= alignmentList.size();
        }
        if (mCnpAlignmentCount != 0) {
            throw new RuntimeException("Unsynchronized alignment count");
        }
    }

    private void writeAlignments(String cnpId, String sampleId, List<Alignment> alignmentList, File inputFile)
        throws IOException {
        File outputDir = new File(".");
        if (mOutputDirectory != null) {
            outputDir = new File(mOutputDirectory);
        }
        String cnpSample = cnpId;
        if (sampleId != null) {
            cnpSample = cnpSample + "_" + sampleId;
        }
        File cnpSampleDir = new File(outputDir, cnpSample);
        if (!cnpSampleDir.exists()) {
            if (!cnpSampleDir.mkdir()) {
                throw new RuntimeException("Failed to create directory " + cnpSampleDir);
            }
        }
        String fileName = inputFile.getName();
        File alignmentFile = new File(cnpSampleDir, fileName);
        PrintWriter writer = new PrintWriter(new FileWriter(alignmentFile, true));
        for (Alignment alignment : alignmentList) {
            writer.println(alignment.arachneFormat());
        }
        writer.flush();
        writer.close();
    }

    private GatherAlignments() {
    }

    private static class CnpRegion {

        private CnpRegion(String cnpId, int sequenceId, int start, int end) {
            mCnpId = cnpId;
            mSequenceId = sequenceId;
            mStart = start;
            mEnd = end;
        }

        public String getCnpId() { return mCnpId; };
        public int getSequenceId() { return mSequenceId; };
        public int getStart() { return mStart; };
        public int getEnd() { return mEnd; };

        private String mCnpId;
        private int mSequenceId;
        private int mStart;
        private int mEnd;
    }

    private boolean mDebug = false;
    private boolean mVerbose = false;

    private boolean mReturnBestHits = false;
    private String mCnpListPath = null;
    private String mSampleId = null;
    private String mInputFileListPath = null;
    private String mOutputDirectory = null;
    private int mCnpRegionPadding = 0;

    private Alignment mPendingAlignment = null;
    private int mCnpAlignmentCount = 0;
    private int mCnpAlignmentLimit = 1000000;
    private Map<String, List<Alignment>> mCnpAlignmentMap = new LinkedHashMap<String, List<Alignment>>();
}



