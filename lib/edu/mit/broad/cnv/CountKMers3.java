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

import java.io.*;
import java.util.*;


/**
 * Tool for counting unique kmers.
 */
public class CountKMers3
{
    private static final int NONUNIQUE_MARKER = -1;
    private static boolean mUseOldFormat = false;

    private String mAction = null;
    private static int mK = 0;
    private int mBatchSize = 0;
    private List<File> mInputFiles = null;
    private File mInputDirectory = null;
    private File mOutputDirectory = null;
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
        new CountKMers3().run(args);
    }

    private void usage() {
        System.out.println("Usage: CountKMers ...");
        System.out.println("  -action <action>");
        System.out.println("  -genome <fasta-file>");
        System.out.println("  -k <k>");
        System.out.println("  -batchSize <n>");
        System.out.println("  -inputDir <directory>");
        System.out.println("  -outputDir <directory>");
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
            } else if (arg.equals("-k") && argsleft > 1) {
                argpos++;
                mK = Integer.parseInt(args[argpos++]);
            } else if (arg.equals("-batchSize") && argsleft > 1) {
                argpos++;
                mBatchSize = Integer.parseInt(args[argpos++]);
            } else if (arg.equals("-inputDir") && argsleft > 1) {
                argpos++;
                mInputDirectory = new File(args[argpos++]);
            } else if (arg.equals("-outputDir") && argsleft > 1) {
                argpos++;
                mOutputDirectory = new File(args[argpos++]);
            } else if (arg.equals("-oldFormat")) {
                argpos++;
                mUseOldFormat = true;
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
            mapKMers();
        } else if (mAction.equals("mapGaps")) {
            mapGaps();
        }
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
        int mapSize = ((mBaseIndex >> 2) & 0x3FFFFFFF) + 1;
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
        if (mUseOldFormat) {
            return readKMerPositionOldFormat(stream);
        }
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
        return new KMerPositionN(encoding, baseIndex);
    }

    private KMerPosition readKMerPositionOldFormat(InputStream stream)
        throws IOException {
        byte[] buffer = mIOBuffer;
        int length = (mK >= 32 ? 20 : 12);
        int count = readFully(stream, buffer, 0, length);
        if (count <= 0) {
            return null;
        } else if (count != length) {
            throw new RuntimeException("Unexpected end of file");
        }
        long encoding = (((long)(buffer[0] & 0xFF)) |
                         ((long)(buffer[1] & 0xFF)) << 8 |
                         ((long)(buffer[2] & 0xFF)) << 16 |
                         ((long)(buffer[3] & 0xFF)) << 24 |
                         ((long)(buffer[4] & 0xFF)) << 32 |
                         ((long)(buffer[5] & 0xFF)) << 40 |
                         ((long)(buffer[6] & 0xFF)) << 48 |
                         ((long)(buffer[7] & 0xFF)) << 56);
        int baseIndex = ((buffer[length-4] & 0xFF) |
                         (buffer[length-3] & 0xFF) << 8 |
                         (buffer[length-2] & 0xFF) << 16 |
                         (buffer[length-1] & 0xFF) << 24);
        if (length == 12) {
            return new KMerPosition1(encoding, baseIndex);
        } else {
            long encoding2 = (((long)(buffer[8] & 0xFF)) |
                              ((long)(buffer[9] & 0xFF)) << 8 |
                              ((long)(buffer[10] & 0xFF)) << 16 |
                              ((long)(buffer[11] & 0xFF)) << 24 |
                              ((long)(buffer[12] & 0xFF)) << 32 |
                              ((long)(buffer[13] & 0xFF)) << 40 |
                              ((long)(buffer[14] & 0xFF)) << 48 |
                              ((long)(buffer[15] & 0xFF)) << 56);
            return new KMerPosition2(encoding, encoding2, baseIndex);
        }
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

    private void writeKMerPosition(OutputStream stream, KMerPosition kmer)
        throws IOException {
        if (mUseOldFormat) {
            writeKMerPositionOldFormat(stream, kmer);
            return;
        }
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

    private void writeKMerPositionOldFormat(OutputStream stream, KMerPosition kmer)
        throws IOException {
        byte[] buffer = mIOBuffer;
        long encoding1 = kmer.getKMerEncoding1();
        long encoding2 = kmer.getKMerEncoding2();
        int baseIndex = kmer.getBaseIndex();
        int offset = 0;
        buffer[offset++] = (byte) ((encoding1) & 0xFF);
        buffer[offset++] = (byte) ((encoding1 >> 8) & 0xFF);
        buffer[offset++] = (byte) ((encoding1 >> 16) & 0xFF);
        buffer[offset++] = (byte) ((encoding1 >> 24) & 0xFF);
        buffer[offset++] = (byte) ((encoding1 >> 32) & 0xFF);
        buffer[offset++] = (byte) ((encoding1 >> 40) & 0xFF);
        buffer[offset++] = (byte) ((encoding1 >> 48) & 0xFF);
        buffer[offset++] = (byte) ((encoding1 >> 56) & 0xFF);
        if (mK >= 32) {
            buffer[offset++] = (byte) ((encoding2) & 0xFF);
            buffer[offset++] = (byte) ((encoding2 >> 8) & 0xFF);
            buffer[offset++] = (byte) ((encoding2 >> 16) & 0xFF);
            buffer[offset++] = (byte) ((encoding2 >> 24) & 0xFF);
            buffer[offset++] = (byte) ((encoding2 >> 32) & 0xFF);
            buffer[offset++] = (byte) ((encoding2 >> 40) & 0xFF);
            buffer[offset++] = (byte) ((encoding2 >> 48) & 0xFF);
            buffer[offset++] = (byte) ((encoding2 >> 56) & 0xFF);
        }
        buffer[offset++] = (byte) ((baseIndex) & 0xFF);
        buffer[offset++] = (byte) ((baseIndex >> 8) & 0xFF);
        buffer[offset++] = (byte) ((baseIndex >> 16) & 0xFF);
        buffer[offset++] = (byte) ((baseIndex >> 24) & 0xFF);
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

    private void createMapFile(int mapSize,
                               File kmerFile,
                               List<StringKMerPosition> exceptionList,
                               File priorMapFile,
                               File mapFile)
        throws IOException {
        byte[] map = null;
        long uniquePriorCount = 0;
        if (priorMapFile.exists()) {
            map = readMapFile(priorMapFile);
            if (map.length != mapSize) {
                throw new RuntimeException("Prior map is wrong size");
            }
            // Clear the new bits from prior map.
            // Also count the prior unique positions while we are at it.
            // Note that this is a count of positions, not kmers.
            for (int i = 0; i < mapSize; i++) {
                int cumBits = map[i] & 0x55;
                uniquePriorCount += Integer.bitCount(cumBits);
                map[i] = (byte) cumBits;
            }
        } else {
            map = new byte[mapSize];
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

        long testCum = 0;
        for (int i = 0; i < map.length; i++) {
            testCum += Integer.bitCount(map[i] & 0x55);
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
        int mod = baseIndex & 0x3;
        int offset = (baseIndex >> 2) & 0x3FFFFFFF;
        if (((map[offset] >> (2*mod)) & 0x3) != 0) {
            throw new RuntimeException("Map entry already set: " + baseIndex);
        }
        map[offset] |= (0x3 << (2*mod));
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
        int byteOffset = (baseIndex >> 2) & 0x3FFFFFFF;
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
        int mod = baseIndex & 0x3;
        return (((mPriorMapValue >> (2*mod)) & 1) != 0);
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

    private static void dbg(String text) {
        System.out.println("#DBG: " + text);
    }

    private static KMerPosition encodeKMer(char[] kmerChars, int baseIndex) {
        if (mUseOldFormat) {
            return encodeKMerOldFormat(kmerChars, baseIndex);
        }
        char[] encoding = encodeKMerChars(kmerChars);
        if (encoding == null) {
            return null;
        }
        char[] reverseEncoding = encodeKMerChars(reverseComplement(kmerChars));
        if (compareEncodings(encoding, reverseEncoding) <= 0) {
            return new KMerPositionN(encoding, baseIndex);
        } else {
            KMerPositionN kmp = new KMerPositionN(reverseEncoding, baseIndex);
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

    private static KMerPosition encodeKMerOldFormat(char[] kmerChars, int baseIndex) {
        if (kmerChars == null) {
            return null;
        }
        int length = kmerChars.length;
        if (length <= 31) {
            long bits = encodeKMerBitsLong(kmerChars, 0, length);
            if (bits == -1) {
                return null;
            }
            return new KMerPosition1(bits, baseIndex);
        } else if (length <= 62) {
            long bits1 = encodeKMerBitsLong(kmerChars, 0, 31);
            long bits2 = encodeKMerBitsLong(kmerChars, 31, length - 31);
            if (bits1 == -1 || bits2 == -1) {
                return null;
            }
            return new KMerPosition2(bits1, bits2, baseIndex);
        } else {
            return null;
        }
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

    private static long encodeKMerBitsLong(char[] kmerChars, int offset, int length) {
        long bits = 0;
        for (int i = 0; i < length; i++) {
            char base = kmerChars[offset + i];
            int baseBits = "ACGT".indexOf(base);
            if (baseBits < 0) {
                return -1;
            }
            bits |= ((long)baseBits) << (2*(length-i-1));
        }
        return bits;
    }

    private static String decodeKMer1(long bits) {
        int length = mK;
        char[] buffer = new char[length];
        decodeKMerBits(bits, buffer, 0, length);
        return new String(buffer);
    }

    private static String decodeKMer2(long bits1, long bits2) {
        int length = mK;
        char[] buffer = new char[length];
        decodeKMerBits(bits1, buffer, 0, 31);
        decodeKMerBits(bits2, buffer, 31, length-31);
        return new String(buffer);
    }

    private static String decodeKMerN(char[] encoding, boolean reverse) {
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

        KMerPosition(int baseIndex) {
            mBaseIndex = baseIndex;
        }

        public String getKMer() {
            return null;
        }

        public long getKMerEncoding1() {
            return -1;
        }

        public long getKMerEncoding2() {
            return -1;
        }

        public final int getBaseIndex() {
            return mBaseIndex;
        }

        public final void setBaseIndex(int baseIndex) {
            mBaseIndex = baseIndex;
        }

        public char[] getKMerEncoding() {
            return null;
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
                   " " + Integer.toHexString(mBaseIndex));
        }
    }

    static class KMerPosition1
        extends KMerPosition {

        private long mKMerEncoding1;

        KMerPosition1(long kmer, int baseIndex) {
            super(baseIndex);
            mKMerEncoding1 = kmer;
        }

        public String getKMer() {
            return decodeKMer1(getKMerEncoding1());
        }

        public final long getKMerEncoding1() {
            return mKMerEncoding1;
        }

        public int compareTo(KMerPosition kmp) {
            int result = Long.signum(getKMerEncoding1() - kmp.getKMerEncoding1());
            if (result == 0) {
                result = Long.signum(getKMerEncoding2() - kmp.getKMerEncoding2());
            }
            return result;
        }
    }

    static class KMerPosition2
        extends KMerPosition1 {

        private long mKMerEncoding2;

        KMerPosition2(long encoding1, long encoding2, int baseIndex) {
            super(encoding1, baseIndex);
            mKMerEncoding2 = encoding2;
        }

        public String getKMer() {
            return decodeKMer2(getKMerEncoding1(), getKMerEncoding2());
        }

        public final long getKMerEncoding2() {
            return mKMerEncoding2;
        }
    }

    static class KMerPositionN
        extends KMerPosition {

        private boolean mReversed;
        private char[] mKMerEncoding;

        KMerPositionN(char[] encoding, int baseIndex) {
            super(baseIndex);
            mReversed = false;
            mKMerEncoding = encoding;
        }

        public boolean getIsReversed() {
            return mReversed;
        }

        public void setIsReversed(boolean value) {
            mReversed = value;
        }

        public String getKMer() {
            return decodeKMerN(mKMerEncoding, mReversed);
        }

        public final char[] getKMerEncoding() {
            return mKMerEncoding;
        }

        public String format() {
            return(getKMer() +
                   " " + formatEncoding(getKMerEncoding()) +
                   " " + (mReversed ? 'R' : 'F') +
                   " " + Integer.toHexString(getBaseIndex()));
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
