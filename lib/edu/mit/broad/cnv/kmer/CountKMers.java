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


import edu.mit.broad.cnv.util.SequenceIterator;

import java.io.*;
import java.util.*;


/**
 * Tool for counting unique kmers.
 */
public class CountKMers
{
    private static final int NONUNIQUE_MARKER = -1;

    private String mAction = null;
    private static int mK = 0;
    private int mMinimumK = 0;
    private int mMaximumK = 0;
    private int mBatchSize = 0;
    private List<File> mInputFiles = null;
    private File mSearchFile = null;
    private String mSequenceName = null;
    private File mInputDirectory = null;
    private File mOutputDirectory = null;
    private boolean mRunDistributed = false;
    private int mDistributedWorkerCount = 0;
    private boolean mVerbose = false;
    private boolean mDebug = false;

    private List<String> mSequenceList = null;
    private List<Integer> mSequenceOffsetList = null;
    private List<File> mSpillFileList = null;
    private double mSpillFactor = 0.9;

    private long mKMerCount = 0;
    private long mUniquePriorCount = 0;
    private long mUniqueNewCount = 0;
    private long mPriorMapUniqueCount = 0;

    private InputStream mPriorMapStream = null;
    private int mPriorMapPosition = -1;
    private int mPriorMapValue = 0;
    private int mInputFileIndex = 0;
    private LineNumberReader mCurrentReader = null;
    private String mNextSequence = null;
    private char[] mKMerBuffer = null;
    private int mKMerBufferedCount = 0;
    private String mLineBuffer = null;
    private int mLineBufferIndex = 0;
    private int mBaseIndex = -1;
    private byte[] mIOBuffer = null;

    /* Design
       Inputs:
       - One or more fasta files to search (currently one).
       - Output directory for the result files.
       - Optionally an input k-1-mer file (output from previous pass).
       Outputs:
       - Unique kmer file: <kmer> <chr> <pos> (sorted by kmer)
         This is unique globally or unique wrt unique (K-1) mers (i.e. K unique, K-1 not).
       - Per chromosome bit map: pos (implicit) new-bit cum-bit
         New-bit is 1 if Kmer starting at pos is unique but (K-1)-mer is not.
         Cum-bit is 1 if Kmer starting at pos is unique for some L <= K.
       - Statistics
       Plan:
       - Reducing memory footprint is crucial.
       - Sequential pass over the input sequences to generate kmers.
       - BatchSize kmers are cached in memory, then sorted and uniqified.
       - As batch array fills, batches are spilled to disk.
       - Batches are reloaded from disk and merged (N-finger algorithm)
       - and streamed to a merge file.
       - Merge file is read from disk and processed as final results.
     */

    public static void main(String[] args)
        throws Exception {
        new CountKMers().run(args);
    }

    private void usage() {
        System.out.println("Usage: CountKMers ...");
        System.out.println("  -action <action>");
        System.out.println("  -genome <fasta-file>");
        System.out.println("  -chromosome <name>");
        System.out.println("  -k <k>");
        System.out.println("  -minK <k>");
        System.out.println("  -maxK <k>");
        System.out.println("  -batchSize <n>");
        System.out.println("  -inputDir <directory>");
        System.out.println("  -outputDir <directory>");
        System.out.println("  -distributed");
        System.out.println("  -workers <n>");
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
            } else if (arg.equals("-genome") && argsleft > 1) {
                argpos++;
                if (mInputFiles == null) {
                    mInputFiles = new ArrayList<File>();
                }
                mInputFiles.add(new File(args[argpos++]));
            } else if (arg.equals("-chromosome") && argsleft > 1) {
                argpos++;
                mSequenceName = args[argpos++];
            } else if (arg.equals("-k") && argsleft > 1) {
                argpos++;
                mK = Integer.parseInt(args[argpos++]);
            } else if (arg.equals("-minK") && argsleft > 1) {
                argpos++;
                mMinimumK = Integer.parseInt(args[argpos++]);
            } else if (arg.equals("-maxK") && argsleft > 1) {
                argpos++;
                mMaximumK = Integer.parseInt(args[argpos++]);
            } else if (arg.equals("-batchSize") && argsleft > 1) {
                argpos++;
                mBatchSize = Integer.parseInt(args[argpos++]);
            } else if (arg.equals("-inputDir") && argsleft > 1) {
                argpos++;
                mInputDirectory = new File(args[argpos++]);
            } else if (arg.equals("-outputDir") && argsleft > 1) {
                argpos++;
                mOutputDirectory = new File(args[argpos++]);
            } else if (arg.equals("-searchFile") && argsleft > 1) {
                argpos++;
                mSearchFile = new File(args[argpos++]);
            } else if (arg.equals("-distributed")) {
                argpos++;
                mRunDistributed = true;
            } else if (arg.equals("-workers") && argsleft > 1) {
                argpos++;
                mDistributedWorkerCount = Integer.parseInt(args[argpos++]);
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
        if (mAction == null || mAction.equals("mapKMers")) {
            if (mRunDistributed) {
                mapKMersDistributed();
            } else {
                mapKMers();
            }
        } else if (mAction.equals("mapGaps")) {
            mapGaps();
        } else if (mAction.equals("rollUp")) {
            rollUp();
        } else if (mAction.equals("search")) {
            search();
        }
    }

    private void search()
        throws IOException {
        char[][] searchStrings = loadSearchFile(mSearchFile);
        while (true) {
            String seqName = getNextSequence();
            if (seqName == null) {
                break;
            }
            int position = 0;
            log("Scanning " + seqName + " ...");
            while (true) {
                char[] kmerChars = getNextKMer();
                if (kmerChars == null) {
                    break;
                }
                position++;
                for (int i = 0; i < searchStrings.length; i++) {
                    if (Arrays.equals(searchStrings[i], kmerChars)) {
                        String kmer = new String(searchStrings[i]);
                        String strand = ((i % 2) == 0) ? "F" : "R";
                        System.out.println(kmer + "\t" + seqName + "\t" + position + "\t" + strand);
                    }
                }
            }
        }
    }

    private char[][] loadSearchFile(File file)
        throws IOException {
        List<char[]> list = new ArrayList<char[]>();
        LineNumberReader reader = new LineNumberReader(new FileReader(file));
        while (true) {
            String line = reader.readLine();
            if (line == null) {
                reader.close();
                break;
            }
            String text = line.trim();
            if (text.length() == 0 || text.startsWith("#")) {
                continue;
            }
            String[] fields = text.split("\\s+");
            char[] kmer = fields[0].toUpperCase().toCharArray();
            list.add(kmer);
            list.add(reverseComplement(kmer));
        }
        return list.toArray(new char[0][0]);
    }

    // Can be used to scan genome for sequence names/lengths.
    private void scanKMers()
        throws IOException {
        mSequenceList = new ArrayList<String>();
        mSequenceOffsetList = new ArrayList<Integer>();
        File priorMapFile =
            new File(mOutputDirectory, "unique_" + (mK-1) + "_mers_map.bin");
        openPriorMap(priorMapFile);
        while (true) {
            String seqName = getNextSequence();
            if (seqName == null) {
                break;
            }
            mSequenceList.add(seqName);
            mSequenceOffsetList.add(mBaseIndex+1);
            log("Scanning " + seqName + " ...");
            while (true) {
                char[] kmerChars = getNextKMer();
                if (kmerChars == null) {
                    break;
                }
                mKMerCount++;
                if (isUniqueInPriorMap(mBaseIndex)) {
                    continue;
                }
            }
        }
        closePriorMap();
    }

    private void mapGaps()
        throws IOException {
        while (true) {
            String seqName = getNextSequence();
            if (seqName == null) {
                break;
            }
            int pos = 0;
            int gapStart = 0;
            while (true) {
                char base = getNextBase();
                if (base == 0) {
                    break;
                }
                pos++;
                if (base == 'N') {
                    if (gapStart == 0) {
                        gapStart = pos;
                    }
                } else {
                    if (gapStart > 0) {
                        System.out.println(seqName + "\t" + gapStart + "\t" + (pos-1));
                        gapStart = 0;
                    }
                }
            }
            if (gapStart > 0) {
                System.out.println(seqName + "\t" + gapStart + "\t" + (pos-1));
                gapStart = 0;
            }
        }
    }

    private void rollUp()
        throws IOException {
        // Roll up based on the middle of the reads.
        File[] mapFiles = getAllMapFiles();
        if (mapFiles.length > 127) {
            throw new RuntimeException("K to large for byte sized counts");
        }
        SequenceIterator seqIterator = new SequenceIterator(mInputFiles);
        while (true) {
            String seqName = seqIterator.getNextSequence();
            if (seqName == null) {
                break;
            }
            if (mSequenceName != null && !mSequenceName.equals(seqName)) {
                continue;
            }
            log("Rolling up sequence " + seqName + " ...");
            int seqBaseIndex = seqIterator.getBaseIndex() + 1;
            char[] seqChars = loadSequence(seqIterator);
            int seqLength = seqChars.length;
            int seqMapOffset = (seqBaseIndex >> 3) & 0x1FFFFFFF;
            int seqMapModulus = (seqBaseIndex & 0x7);
            int seqMapLength = (seqMapModulus + seqLength + 7)/8;
            // log("  seqLength = " + seqLength);
            // log("  baseIndex = " + Integer.toHexString(seqBaseIndex)
            //     + " (" + (((long)seqBaseIndex) & 0xFFFFFFFFL) + ")");
            // log("  seqMapOffset = " + seqMapOffset);
            // log("  seqMapLength = " + seqMapLength);
            byte[] counts = new byte[seqLength];
            for (int pos = 1; pos <= seqLength; pos++) {
                if (seqChars[pos-1] == 'N') {
                    counts[pos-1] = -1;
                }
            }
            for (int k = 1; k <= mapFiles.length; k++) {
                if (mapFiles[k-1] == null) {
                    continue;
                }
                log("Processing map file " + mapFiles[k-1] + " ...");
                byte[] kmerMap = readMapFileRegion(mapFiles[k-1], seqMapOffset, seqMapLength);
                for (int pos = 1; pos <= seqLength; pos++) {
                    if (counts[pos-1] != 0) {
                        continue;
                    } else if (isNearContigBoundary(pos, seqChars, k)) {
                        counts[pos-1] = -1;
                    } else {
                        int baseOffset = pos - (k+1)/2;
                        int mapIndex = seqMapModulus + baseOffset;
                        if (isUniqueInMap(kmerMap, mapIndex)) {
                            counts[pos-1] = (byte) k;
                        }
                    }
                }
            }
            File outputFile =
                new File(mOutputDirectory, "rollup_" + seqName + ".bin");
            writeRollUpFile(outputFile, counts);
        }
    }

    private boolean isNearContigBoundary(int pos, char[] seqChars, int k) {
        int windowStart = pos - (k-1)/2;
        int windowEnd = pos + k/2;
        if (windowStart < 1 || windowEnd > seqChars.length) {
            return true;
        }
        for (int i = windowStart-1; i < windowEnd; i++) {
            if (seqChars[i] == 'N') {
                return true;
            }
        }
        return false;
    }

    private void writeRollUpFile(File file, byte[] counts)
        throws IOException {
        FileOutputStream stream = new FileOutputStream(file);
        stream.write(counts);
        stream.flush();
        stream.close();
        if (mDebug) {
            PrintWriter writer = new PrintWriter(file + ".dbg");
            for (int i = 0; i < counts.length; i++) {
                writer.println(counts[i]);
            }
            writer.flush();
            writer.close();
        }
    }

    /**
     * Returns an array of files, indexed by K,
     * where the array index = K-1 (i.e. K=1 is the first file).
     * If there is no file for index K, then the array element is null.
     */
    private File[] getAllMapFiles() {
        int maxK = mMaximumK;
        if (maxK == 0) {
            // Safe upper bound
            maxK = 1000;
        }
        List<File> fileList = new ArrayList<File>();
        for (int k = 1; k <= maxK; k++) {
            if (mMinimumK > 0 && k < mMinimumK) {
                continue;
            }
            File mapFile =
                new File(mInputDirectory, "unique_" + k + "_mers_map.bin");
            if (mapFile.exists()) {
                while (fileList.size() < k-1) {
                    fileList.add(null);
                }
                fileList.add(mapFile);
            } else {
                if (mMaximumK == 0 && !fileList.isEmpty()) {
                    break;
                }
            }
        }
        File[] result = new File[fileList.size()];
        result = fileList.toArray(result);
        if (mDebug) {
            for (int i = 0; i < result.length; i++) {
                debug("mapFiles[k=" + (i+1) + "] = " + result[i]);
            }
        }
        return result;
    }

    private char[] loadSequence(SequenceIterator seqIterator)
        throws IOException {
        StringBuilder builder = new StringBuilder();
        while (true) {
            char ch = seqIterator.getNextBase();
            if (ch == 0) {
                break;
            }
            builder.append(ch);
        }
        char[] result = new char[builder.length()];
        builder.getChars(0, builder.length(), result, 0);
        return result;
    }

    private void mapKMersDistributed()
        throws Exception {
        DistributedKMerCounter algorithm = new DistributedKMerCounter();
        algorithm.setDebug(mDebug);
        algorithm.setVerbose(mVerbose);
        algorithm.setInputFiles(mInputFiles);
        algorithm.setK(mK);
        algorithm.setMaximumWorkerCount(mDistributedWorkerCount);
        // algorithm.setLsfQueue(mLsfQueue);
        // algorithm.setLsfLogDirectory(mLsfLogDirectory);
        // algorithm.setEnableGcLogging(mEnableGcLogging);
        algorithm.run();
    }

    private void mapKMers()
        throws IOException {

        File textKMerFile =
            new File(mOutputDirectory, "unique_" + mK + "_mers.txt");
        File binaryKMerFile =
            new File(mOutputDirectory, "unique_" + mK + "_mers.bin");
        File exceptionFile =
            new File(mOutputDirectory, "unique_" + mK + "_mers.extra");
        File mapFile =
            new File(mOutputDirectory, "unique_" + mK + "_mers_map.bin");
        File priorMapFile =
            new File(mOutputDirectory, "unique_" + (mK-1) + "_mers_map.bin");
        File statsFile =
            new File(mOutputDirectory, "unique_" + mK + "_mers_stats.txt");

        if (mBatchSize == 0) {
            throw new RuntimeException("Batch size not specified");
        }

        int kmerCount = 0;
        int batchSize = mBatchSize;
        KMerPosition[] kmerArray = new KMerPosition[batchSize];
        List<StringKMerPosition> exceptionList = new ArrayList<StringKMerPosition>();
        mSequenceList = new ArrayList<String>();
        mSequenceOffsetList = new ArrayList<Integer>();
        mIOBuffer = new byte[Math.max(20,4 + 2*((mK + 7)/8))];

        openPriorMap(priorMapFile);

        while (true) {
            String seqName = getNextSequence();
            if (seqName == null) {
                break;
            }
            mSequenceList.add(seqName);
            mSequenceOffsetList.add(mBaseIndex+1);
            log("Processing " + seqName + " ...");
            while (true) {
                char[] kmerChars = getNextKMer();
                if (kmerChars == null) {
                    break;
                }
                mKMerCount++;
                int baseIndex = mBaseIndex;
                if (isUniqueInPriorMap(baseIndex)) {
                    mUniquePriorCount++;
                    continue;
                }
                KMerPosition kmp = encodeKMer(kmerChars, baseIndex);
                if (kmp == null) {
                    // Note: We currently do not handle the reverse
                    // complement of exception characters correctly.
                    // For hg18, however, this doesn't matter as
                    // none of the kmers containing non-ACGT characters
                    // are present on the reverse strand.
                    String kmer = new String(kmerChars);
                    exceptionList.add(new StringKMerPosition(kmer, baseIndex));
                    continue;
                }
                kmerArray[kmerCount++] = kmp;
                if (kmerCount == batchSize) {
                    kmerCount = compactKMers(kmerArray, kmerCount);
                    if (kmerCount > mSpillFactor * batchSize) {
                        spillKMers(kmerArray, kmerCount);
                        kmerCount = 0;
                    }
                }
            }
        }
        if (kmerCount > 0) {
            kmerCount = compactKMers(kmerArray, kmerCount);
            if (mSpillFileList != null) {
                spillKMers(kmerArray, kmerCount);
                kmerCount = 0;
            }
        }

        closePriorMap();

        // Write out the exception kmers (text file).
        compactKMers(exceptionList);
        writeExceptionFile(exceptionList, exceptionFile);

        // Write out the binary file of unique encoded kmers.
        if (mSpillFileList == null) {
            kmerCount = removeNonUnique(kmerArray, kmerCount);
            writeKMerBinaryFile(kmerArray, kmerCount, binaryKMerFile);
            mUniqueNewCount = kmerCount;
        } else {
            mUniqueNewCount = mergeSpillFiles(mSpillFileList, binaryKMerFile);
        }
        mUniqueNewCount += countUniqueKMers(exceptionList);

        // Write out the text file of (all) unique kmers.
        writeKMerTextFile(binaryKMerFile, exceptionList, textKMerFile);

        // Create map file from prior map plus the new unique kmers.
        long mapSize = (mBaseIndex + 1) & 0xFFFFFFFFL;
        createMapFile(mapSize, binaryKMerFile, exceptionList, priorMapFile, mapFile);

        // Write summary statistics file.
        writeSummaryStatistics(statsFile);
    }

    private int compactKMers(KMerPosition[] kmerArray, int kmerCount) {
        if (kmerCount == 0) {
            return 0;
        }
        log("Compacting " + kmerCount + " kmers at index " +
            Integer.toHexString(mBaseIndex) + " ...");
        Arrays.sort(kmerArray, 0, kmerCount);
        int newCount = 1;
        KMerPosition current = kmerArray[0];
        for (int i = 1; i < kmerCount; i++) {
            KMerPosition kmp = kmerArray[i];
            if (current.compareTo(kmp) == 0) {
                current.setBaseIndex(NONUNIQUE_MARKER);
            } else {
                kmerArray[newCount++] = kmp;
                current = kmp;
            }
        }
        log("Compaction finished, new count is " + newCount);
        return newCount;
    }

    private int compactKMers(StringKMerPosition[] kmerArray, int kmerCount) {
        if (kmerCount == 0) {
            return 0;
        }
        log("Compacting " + kmerCount + " string kmers ...");
        Arrays.sort(kmerArray, 0, kmerCount);
        int newCount = 1;
        String kmerString = kmerArray[0].getKMer();
        for (int i = 1; i < kmerCount; i++) {
            StringKMerPosition kmp = kmerArray[i];
            String ks = kmp.getKMer();
            if (ks.equals(kmerString)) {
                kmerArray[newCount-1].setBaseIndex(NONUNIQUE_MARKER);
            } else {
                kmerArray[newCount++] = kmp;
                kmerString = ks;
            }
        }
        log("Compaction finished, new count is " + newCount);
        return newCount;
    }

    private void compactKMers(List<StringKMerPosition> kmerList) {
        int kmerCount = kmerList.size();
        if (kmerCount <= 1) {
            return;
        }
        StringKMerPosition[] kmerArray =
            kmerList.toArray(new StringKMerPosition[kmerCount]);
        kmerCount = compactKMers(kmerArray, kmerCount);
        kmerList.clear();
        for (int i = 0; i < kmerCount; i++) {
            kmerList.add(kmerArray[i]);
        }
    }

    private int removeNonUnique(KMerPosition[] kmerArray, int kmerCount) {
        int uniqueCount = 0;
        for (int i = 0; i < kmerCount; i++) {
            KMerPosition kmp = kmerArray[i];
            if (kmp.getBaseIndex() != NONUNIQUE_MARKER) {
                kmerArray[uniqueCount++] = kmp;
            }
        }
        return uniqueCount;
    }

    private int countUniqueKMers(List<StringKMerPosition> kmerList) {
        int uniqueCount = 0;
        for (StringKMerPosition kmp : kmerList) {
            if (kmp.getBaseIndex() != NONUNIQUE_MARKER) {
                uniqueCount++;
            }
        }
        return uniqueCount;
    }

    private void spillKMers(KMerPosition[] kmerArray, int kmerCount)
        throws IOException {
        if (mSpillFileList == null) {
            mSpillFileList = new ArrayList<File>();
        }
        int fileNumber = mSpillFileList.size() + 1;
        log("Spilling " + kmerCount + " kmers to file " + fileNumber + " ...");
        File spillFile = new File(mOutputDirectory,
                                  "spill_" + mK + "_" + fileNumber + ".tmp");
        mSpillFileList.add(spillFile);
        writeKMerBinaryFile(kmerArray, kmerCount, spillFile);
        log("Spill file written");
    }

    private void writeKMerBinaryFile(KMerPosition[] kmerArray,
                                     int kmerCount,
                                     File outputFile)
        throws IOException {
        OutputStream outputStream =
            new BufferedOutputStream(new FileOutputStream(outputFile));
        for (int i = 0; i < kmerCount; i++) {
            KMerPosition kmp = kmerArray[i];
            writeKMerPosition(outputStream, kmerArray[i]);
        }
        outputStream.flush();
        outputStream.close();
    }

    private void writeExceptionFile(List<StringKMerPosition> kmerList,
                                    File outputFile)
        throws IOException {
        PrintWriter writer =
            new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
        for (StringKMerPosition kmer : kmerList) {
            writeUniqueKMer(kmer, writer);
        }
        writer.flush();
        writer.close();
    }

    private KMerPosition readKMerPosition(InputStream stream)
        throws IOException {
        byte[] buffer = mIOBuffer;
        int encodingLength = (mK + 7)/8;
        int fileLength = 4 + 2*encodingLength;
        int count = readFully(stream, buffer, 0, fileLength);
        if (count <= 0) {
            return null;
        } else if (count != fileLength) {
            throw new RuntimeException("Unexpected end of file");
        }
        char[] encoding = new char[encodingLength];
        int baseIndex = ((buffer[0] & 0xFF) |
                         (buffer[1] & 0xFF) << 8 |
                         (buffer[2] & 0xFF) << 16 |
                         (buffer[3] & 0xFF) << 24);
        for (int i = 0; i < encodingLength; i++) {
            encoding[i] = (char) ((buffer[2*i+4] & 0xFF) |
                                  ((buffer[2*i+5] & 0xFF) << 8));
        }
        return new KMerPosition(encoding, baseIndex);
    }

    private int readFully(InputStream stream, byte[] buffer, int offset, int count)
        throws IOException {
        int readCount = 0;
        while (readCount < count) {
            int read = stream.read(buffer, offset, count-readCount);
            if (read <= 0) {
                break;
            }
            offset += read;
            readCount += read;
        }
        return readCount;
    }

    private void skipBytes(InputStream stream, int count)
        throws IOException {

        long longCount = count;
        long skipCount = 0;
        while (skipCount < longCount) {
            long skipped = stream.skip(longCount - skipCount);
            if (skipped <= 0) {
                throw new RuntimeException("Skip failed");
            }
            skipCount += skipped;
        }
    }

    private void writeKMerPosition(OutputStream stream, KMerPosition kmer)
        throws IOException {
        byte[] buffer = mIOBuffer;
        int baseIndex = kmer.getBaseIndex();
        char[] encoding = kmer.getKMerEncoding();
        int offset = 0;
        buffer[offset++] = (byte) ((baseIndex) & 0xFF);
        buffer[offset++] = (byte) ((baseIndex >> 8) & 0xFF);
        buffer[offset++] = (byte) ((baseIndex >> 16) & 0xFF);
        buffer[offset++] = (byte) ((baseIndex >> 24) & 0xFF);
        for (int i = 0; i < encoding.length; i++) {
            buffer[offset++] = (byte) ((encoding[i]) & 0xFF);
            buffer[offset++] = (byte) ((encoding[i] >> 8) & 0xFF);
        }
        stream.write(buffer, 0, offset);
    }

    private long mergeSpillFiles(List<File> spillFiles, File outputFile)
        throws IOException {

        if (spillFiles == null) {
            return 0;
        }

        log("Merging spill files ...");
        OutputStream outputStream =
            new BufferedOutputStream(new FileOutputStream(outputFile));
        long uniqueCount = 0;
        int fileCount = spillFiles.size();
        InputStream[] inputStreams = new InputStream[fileCount];
        KMerPosition[] kmers = new KMerPosition[fileCount];
        for (int i = 0; i < fileCount; i++) {
            inputStreams[i] =
                new BufferedInputStream(new FileInputStream(spillFiles.get(i)));
        }
        while (true) {
            for (int i = 0; i < fileCount; i++) {
                if (kmers[i] == null && inputStreams[i] != null) {
                    kmers[i] = readKMerPosition(inputStreams[i]);
                    if (kmers[i] == null) {
                        inputStreams[i].close();
                        inputStreams[i] = null;
                    }
                }
            }
            int count = 0;
            KMerPosition kmer = null;
            for (int i = 0; i < fileCount; i++) {
                KMerPosition kmp = kmers[i];
                if (kmp == null) {
                    continue;
                } else if (kmer == null) {
                    kmer = kmp;
                    count = 1;
                } else {
                    int cmp = kmp.compareTo(kmer);
                    if (cmp == 0) {
                        count++;
                    } else if (cmp < 0) {
                        kmer = kmp;
                        count = 1;
                    }
                }
            }
            if (kmer == null) {
                break;
            }
            for (int i = 0; i < fileCount; i++) {
                if (kmers[i] != null && kmer.compareTo(kmers[i]) == 0) {
                    kmers[i] = null;
                }
            }
            if (count == 1 && kmer.getBaseIndex() != NONUNIQUE_MARKER) {
                uniqueCount++;
                writeKMerPosition(outputStream, kmer);
            }
        }
        outputStream.flush();
        outputStream.close();
        for (int i = 0; i < fileCount; i++) {
            // spillFiles.get(i).delete();
        }
        log("Spill files merged, unique count is " + uniqueCount);
        return uniqueCount;
    }

    private void writeKMerTextFile(File inputFile,
                                   List<StringKMerPosition> exceptionList,
                                   File outputFile)
        throws IOException {

        log("Writing kmer file " + outputFile + " ...");
        int exceptionIndex = 0;
        StringKMerPosition excKMer = null;
        Iterator<StringKMerPosition> excIter = null;
        if (!exceptionList.isEmpty()) {
            excIter = exceptionList.iterator();
            excKMer = excIter.next();
        }

        InputStream inputStream =
            new BufferedInputStream(new FileInputStream(inputFile));
        PrintWriter writer =
            new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
        KMerPosition kmer = readKMerPosition(inputStream);
        while (kmer != null || excKMer != null) {
            if (excKMer == null) {
                writeUniqueKMer(kmer, writer);
                kmer = readKMerPosition(inputStream);
            } else if (kmer == null) {
                writeUniqueKMer(excKMer, writer);
                excKMer = excIter.hasNext() ? excIter.next() : null;
            } else if (kmer.getKMer().compareTo(excKMer.getKMer()) < 0) {
                writeUniqueKMer(kmer, writer);
                kmer = readKMerPosition(inputStream);
            } else {
                writeUniqueKMer(excKMer, writer);
                excKMer = excIter.hasNext() ? excIter.next() : null;
            }
        }
        inputStream.close();
        writer.flush();
        writer.close();
        log("Wrote kmer file: " + outputFile);
    }

    private void writeUniqueKMer(KMerPosition kmer, PrintWriter writer) {
        if (kmer.getBaseIndex() != NONUNIQUE_MARKER) {
            writeKMer(kmer.getKMer(), kmer.getBaseIndex(), writer);
        }
    }

    private void writeUniqueKMer(StringKMerPosition kmer, PrintWriter writer) {
        if (kmer.getBaseIndex() != NONUNIQUE_MARKER) {
            writeKMer(kmer.getKMer(), kmer.getBaseIndex(), writer);
        }
    }

    private void writeKMer(String kmer, int baseIndex, PrintWriter writer) {
        String chr = getBaseIndexSequenceName(baseIndex);
        int pos = getBaseIndexCoordinate(baseIndex);
        writer.println(kmer + "\t" + chr + "\t" + pos);
    }

    private void createMapFile(long mapSize,
                               File kmerFile,
                               List<StringKMerPosition> exceptionList,
                               File priorMapFile,
                               File mapFile)
        throws IOException {
        byte[] map = null;
        long uniquePriorCount = 0;
        long byteSize = (mapSize + 7)/8;
        int mapByteSize = (int) byteSize;
        if (mapByteSize != byteSize) {
            throw new RuntimeException("Map too large: " + mapSize);
        }
        if (priorMapFile.exists()) {
            map = readMapFile(priorMapFile);
            if (map.length != mapByteSize) {
                throw new RuntimeException("Prior map is wrong size");
            }
            // Count the prior unique positions
            for (int i = 0; i < mapByteSize; i++) {
                uniquePriorCount += Integer.bitCount(map[i] & 0xFF);
            }
        } else {
            map = new byte[mapByteSize];
        }
        for (StringKMerPosition kmp : exceptionList) {
            addToMap(kmp, map);
        }
        mPriorMapUniqueCount = uniquePriorCount;

        InputStream inputStream =
            new BufferedInputStream(new FileInputStream(kmerFile));
        while (true) {
            KMerPosition kmp = readKMerPosition(inputStream);
            if (kmp == null) {
                inputStream.close();
                break;
            }
            addToMap(kmp, map);
        }

        writeMapFile(map, mapFile);
    }

    private void addToMap(KMerPosition kmp, byte[] map) {
        int baseIndex = kmp.getBaseIndex();
        if (baseIndex != NONUNIQUE_MARKER) {
            addToMap(baseIndex, map);
        }
    }

    private void addToMap(StringKMerPosition kmp, byte[] map) {
        int baseIndex = kmp.getBaseIndex();
        if (baseIndex != NONUNIQUE_MARKER) {
            addToMap(baseIndex, map);
        }
    }

    private void addToMap(int baseIndex, byte[] map) {
        int mod = baseIndex & 0x7;
        int offset = (baseIndex >> 3) & 0x1FFFFFFF;
        if ((map[offset] & (1 << mod)) != 0) {
            throw new RuntimeException("Map entry already set: " + baseIndex);
        }
        map[offset] |= (1 << mod);
    }

    private boolean isUniqueInMap(byte[] map, int baseIndex) {
        int mod = baseIndex & 0x7;
        int offset = (baseIndex >> 3) & 0x1FFFFFFF;
        return ((map[offset] & (1 << mod)) != 0);
    }

    private void writeSummaryStatistics(File outputFile)
        throws IOException {
        PrintWriter writer =
            new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
        long baseCount = (mBaseIndex + 1) & 0xFFFFFFFFL;
        long uniqueCount = mUniquePriorCount + mUniqueNewCount;
        long nonUniqueCount = mKMerCount - uniqueCount;
        writer.println("K: " + mK);
        writer.println("Sequences: " + mSequenceList.size());
        writer.println("Bases: " + baseCount);
        writer.println("KMers: " + mKMerCount);
        writer.println("Prior map count: " + mPriorMapUniqueCount);
        writer.println("Unique prior: " + mUniquePriorCount +
                       " (" + formatPercent(mUniquePriorCount, mKMerCount) + ")");
        writer.println("Unique new: " + mUniqueNewCount +
                       " (" + formatPercent(mUniqueNewCount, mKMerCount) + ")");
        writer.println("Unique cumulative: " + uniqueCount +
                       " (" + formatPercent(uniqueCount, mKMerCount) + ")");
        writer.println("Nonunique: " + nonUniqueCount +
                       " (" + formatPercent(nonUniqueCount, mKMerCount) + ")");
        writer.flush();
        writer.close();
    }

    private String formatPercent(long numerator, long denominator) {
        double fraction = 0.0;
        if (denominator != 0) {
            fraction = numerator / (double) denominator;
        }
        return String.format("%1.1f%%", fraction * 100.0);
    }

    private void openPriorMap(File mapFile)
        throws IOException {
        if (mapFile.exists()) {
            mPriorMapStream = new BufferedInputStream(new FileInputStream(mapFile));
            mPriorMapPosition = -1;
            mPriorMapValue = 0;
        }
    }

    private void closePriorMap()
        throws IOException {
        if (mPriorMapStream != null) {
            mPriorMapStream.close();
        }
        mPriorMapStream = null;
        mPriorMapPosition = -1;
        mPriorMapValue = 0;
    }

    private byte[] readMapFile(File file)
        throws IOException {
        long fileLength = file.length();
        if (fileLength > 1000000000) {
            throw new RuntimeException("Prior map too large: " + file);
        }
        int length = (int) fileLength;
        byte[] map = new byte[length];
        FileInputStream stream = new FileInputStream(file);
        int count = readFully(stream, map, 0, length);
        if (count != length) {
            throw new RuntimeException("Failed to read map: " + file);
        }
        stream.close();
        return map;
    }

    /**
     * Read just a subset of a map file.
     */
    private byte[] readMapFileRegion(File file, int offset, int length)
        throws IOException {
        byte[] map = new byte[length];
        FileInputStream stream = new FileInputStream(file);
        skipBytes(stream, offset);
        int count = readFully(stream, map, 0, length);
        if (count != length) {
            throw new RuntimeException("Failed to read map: " + file);
        }
        stream.close();
        return map;
    }

    private void writeMapFile(byte[] map, File file)
        throws IOException {
        FileOutputStream stream = new FileOutputStream(file);
        stream.write(map);
        stream.flush();
        stream.close();
    }

    private boolean isUniqueInPriorMap(int baseIndex)
        throws IOException {
        if (mPriorMapStream == null) {
            return false;
        }
        int byteOffset = (baseIndex >> 3) & 0x1FFFFFFF;
        if (byteOffset != mPriorMapPosition) {
            int delta = byteOffset - mPriorMapPosition;
            if (delta < 0) {
                throw new RuntimeException("Attempt to seek backwards in prior map");
            }
            if (delta > 1) {
                skipFully(mPriorMapStream, delta-1);
            }
            mPriorMapValue = mPriorMapStream.read();
            if (mPriorMapValue < 0) {
                throw new RuntimeException("Unexpected end of file in prior map");
            }
            mPriorMapPosition += delta;
        }
        int mod = baseIndex & 0x7;
        return (((1 << mod) & mPriorMapValue) != 0);
    }

    private void skipFully(InputStream stream, long amount)
        throws IOException {
        while (amount > 0) {
            long skip = stream.skip(amount);
            if (skip <= 0 || skip > amount) {
                throw new RuntimeException("Skip failed");
            }
            amount -= skip;
        }
    }

    private String getBaseIndexSequenceName(int baseIndex) {
        int sequenceCount = mSequenceList.size();
        for (int i = 0; i < sequenceCount-1; i++) {
            int nextOffset = mSequenceOffsetList.get(i+1);
            if (compareBaseIndex(nextOffset, baseIndex) > 0) {
                return mSequenceList.get(i);
            }
        }
        return mSequenceList.get(sequenceCount-1);
    }

    private int getBaseIndexCoordinate(int baseIndex) {
        Integer sequenceOffset = null;
        for (Integer offset : mSequenceOffsetList) {
            if (compareBaseIndex(offset, baseIndex) > 0) {
                break;
            }
            sequenceOffset = offset;
        }
        if (sequenceOffset == null) {
            return 0;
        }
        int coordinate = baseIndex - sequenceOffset + 1;
        if (coordinate <= 0) {
            dumpSequenceList();
            System.out.println("coordinate: " + coordinate);
            System.out.println("sequenceOffset: " + Integer.toHexString(sequenceOffset));
            System.out.println("baseIndex: " + Integer.toHexString(baseIndex));
            throw new RuntimeException("Internal error: illegal coordinate " +
                                       coordinate + " for base index " + baseIndex);
        }
        return coordinate;
    }

    private void dumpSequenceList() {
        System.out.println("# Sequences:");
        int count = mSequenceList.size();
        for (int i = 0; i < count; i++) {
            String seqName = mSequenceList.get(i);
            int offset = mSequenceOffsetList.get(i);
            System.out.println("# " + seqName +
                               "\t" + offset +
                               "\t" + Integer.toHexString(offset));
        }
    }

    private int compareBaseIndex(int baseIndex1, int baseIndex2) {
        // Implements unsigned comparison, a la compareTo
        if (baseIndex1 < 0 ^ baseIndex2 < 0) {
            return ((baseIndex1 < 0) ? 1 : -1);
        } else {
            return (baseIndex1 - baseIndex2);
        }
    }

    private String getNextSequence()
        throws IOException {

        while (mNextSequence == null) {
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
            }
        }
        String result = mNextSequence;
        mNextSequence = null;
        return result;
    }

    private LineNumberReader getNextReader()
        throws IOException {
        if (mInputFileIndex >= mInputFiles.size()) {
            return null;
        }
        File file = mInputFiles.get(mInputFileIndex++);
        return new LineNumberReader(new FileReader(file));
    }

    private char[] getNextKMer()
        throws IOException {

        if (mKMerBuffer == null) {
            mKMerBuffer = new char[mK];
        }
        System.arraycopy(mKMerBuffer, 1, mKMerBuffer, 0, mKMerBuffer.length - 1);
        if (mKMerBufferedCount > 0) {
            mKMerBufferedCount--;
        }

        while (mKMerBufferedCount < mK) {
            char base = getNextBase();
            if (base == 0) {
                incrementBaseIndex(mKMerBufferedCount);
                mKMerBufferedCount = 0;
                return null;
            } else if (base == 'N') {
                incrementBaseIndex(mKMerBufferedCount+1);
                mKMerBufferedCount = 0;
            } else {
                mKMerBuffer[mKMerBufferedCount++] = base;
            }
        }
        incrementBaseIndex(1);
        return mKMerBuffer;
    }

    private char getNextBase()
        throws IOException {

        if (mLineBuffer == null || mLineBufferIndex >= mLineBuffer.length()) {
            if (mCurrentReader == null) {
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
        return mLineBuffer.charAt(mLineBufferIndex++);
    }

    private void incrementBaseIndex(int amount) {
        if (mBaseIndex < -1 && (mBaseIndex + amount) >= -1) {
            throw new RuntimeException("Base index: 32-bit overflow");
        }
        mBaseIndex += amount;
    }

    private void log(String text) {
        if (mVerbose) {
            System.out.println("# " + new Date() + " " + text);
        }
    }

    private void debug(String text) {
        if (mDebug) {
            System.out.println("# " + new Date() + " " + text);
        }
    }

    private static KMerPosition encodeKMer(char[] kmerChars, int baseIndex) {
        char[] encoding = encodeKMerChars(kmerChars);
        if (encoding == null) {
            return null;
        }
        char[] reverseEncoding = encodeKMerChars(reverseComplement(kmerChars));
        if (compareEncodings(encoding, reverseEncoding) <= 0) {
            return new KMerPosition(encoding, baseIndex);
        } else {
            KMerPosition kmp = new KMerPosition(reverseEncoding, baseIndex);
            kmp.setIsReversed(true);
            return kmp;
        }
    }

    private static char[] encodeKMerChars(char[] kmerChars) {
        if (kmerChars == null) {
            return null;
        }

        int kmerLength = kmerChars.length;
        int encodingLength = (kmerLength + 7) / 8;
        char[] encoding = new char[encodingLength];
        int offset = kmerLength % 8;
        offset = (offset == 0) ? 8 : offset;
        int bits = encodeKMerBits(kmerChars, 0, offset);
        if (bits < 0) {
            return null;
        }
        encoding[0] = (char) bits;
        for (int i = 1; i < encodingLength; i++) {
            bits = encodeKMerBits(kmerChars, offset, 8);
            if (bits < 0) {
                return null;
            }
            encoding[i] = (char) bits;
            offset += 8;
        }
        return encoding;
    }

    private static int compareEncodings(char[] encoding1, char[] encoding2) {
        int length = Math.max(encoding1.length, encoding2.length);
        for (int i = 0; i < length; i++) {
            int result = encoding1[i] - encoding2[i];
            if (result != 0) {
                return result;
            }
        }
        return 0;
    }

    private static int encodeKMerBits(char[] kmerChars, int offset, int length) {
        int bits = 0;
        for (int i = 0; i < length; i++) {
            char base = kmerChars[offset + i];
            int baseBits = "ACGT".indexOf(base);
            if (baseBits < 0) {
                return -1;
            }
            bits |= baseBits << (2*(length-i-1));
        }
        return bits;
    }

    private static String decodeKMer(char[] encoding, boolean reverse) {
        int length = mK;
        char[] buffer = new char[length];
        int offset = length % 8;
        offset = (offset == 0) ? 8 : offset;
        decodeKMerBits(encoding[0], buffer, 0, offset);
        for (int i = 1; i < encoding.length; i++) {
            decodeKMerBits(encoding[i], buffer, offset, 8);
            offset += 8;
        }
        if (reverse) {
            reverseComplementInPlace(buffer);
        }
        return new String(buffer);
    }

    private static void decodeKMerBits(char bits, char[] buffer, int offset, int length) {
        for (int i = 0; i < length; i++) {
            int baseBits = (int) ((bits >> (2*(length-i-1))) & 0x3);
            buffer[offset + i] = "ACGT".charAt(baseBits);
        }
    }

    private static void decodeKMerBits(long bits, char[] buffer, int offset, int length) {
        for (int i = 0; i < length; i++) {
            int baseBits = (int) ((bits >> (2*(length-i-1))) & 0x3);
            buffer[offset + i] = "ACGT".charAt(baseBits);
        }
    }

    private static char[] reverseComplement(char[] buffer) {
        int length = buffer.length;
        char[] result = new char[length];
        System.arraycopy(buffer, 0, result, 0, length);
        reverseComplementInPlace(result);
        return result;
    }

    private static void reverseComplementInPlace(char[] buffer) {
        int length = buffer.length;
        int limit = (length + 1)/2;
        for (int i = 0; i < limit; i++) {
            char ch1 = reverseComplement(buffer[i]);
            char ch2 = reverseComplement(buffer[length-i-1]);
            buffer[i] = ch2;
            buffer[length-i-1] = ch1;
        }
    }

    private static char reverseComplement(char base) {
        switch (base) {
            case 'A':
                return 'T';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            case 'T':
                return 'A';
        }
        return base;
    }

    private static String formatEncoding(char[] encoding) {
        if (encoding == null) {
            return null;
        }
        StringBuilder builder = new StringBuilder();
        builder.append('[');
        for (int i = 0; i < encoding.length; i++) {
            String hex = Integer.toHexString(encoding[i]);
            int length = hex.length();
            while (length < 4) {
                builder.append('0');
                length++;
            }
            builder.append(hex);
        }
        builder.append(']');
        return builder.toString();
    }

    static class KMerPosition
        implements Comparable<KMerPosition> {

        private int mBaseIndex;
        private boolean mReversed;
        private char[] mKMerEncoding;

        KMerPosition(char[] encoding, int baseIndex) {
            mBaseIndex = baseIndex;
            mReversed = false;
            mKMerEncoding = encoding;
        }

        public final String getKMer() {
            return decodeKMer(mKMerEncoding, mReversed);
        }

        public final boolean getIsReversed() {
            return mReversed;
        }

        public final void setIsReversed(boolean value) {
            mReversed = value;
        }

        public final int getBaseIndex() {
            return mBaseIndex;
        }

        public final void setBaseIndex(int baseIndex) {
            mBaseIndex = baseIndex;
        }

        public final char[] getKMerEncoding() {
            return mKMerEncoding;
        }

        public int compareTo(KMerPosition kmp) {
            return compareEncodings(getKMerEncoding(), kmp.getKMerEncoding());
        }

        public boolean equals(Object object) {
            if (!(object instanceof KMerPosition)) {
                return false;
            }
            KMerPosition kmp = (KMerPosition) object;
            return (getBaseIndex() == kmp.getBaseIndex() &&
                    this.compareTo(kmp) == 0);
        }

        public String format() {
            return(getKMer() +
                   " " + formatEncoding(getKMerEncoding()) +
                   " " + (mReversed ? 'R' : 'F') +
                   " " + Integer.toHexString(mBaseIndex));
        }
    }

    static class StringKMerPosition
        implements Comparable<StringKMerPosition> {

        private String mKMerString = null;
        private int mBaseIndex;

        StringKMerPosition(String kmer, int baseIndex) {
            mKMerString = kmer;
            mBaseIndex = baseIndex;
        }

        public final String getKMer() {
            return mKMerString;
        }

        public final int getBaseIndex() {
            return mBaseIndex;
        }

        public final void setBaseIndex(int baseIndex) {
            mBaseIndex = baseIndex;
        }

        public int compareTo(StringKMerPosition kmp) {
            return mKMerString.compareTo(kmp.mKMerString);
        }

        public boolean equals(Object object) {
            if (!(object instanceof StringKMerPosition)) {
                return false;
            }
            StringKMerPosition kmp = (StringKMerPosition) object;
            return (mBaseIndex == kmp.mBaseIndex &&
                    mKMerString.equals(kmp.mKMerString));
        }
    }
}
