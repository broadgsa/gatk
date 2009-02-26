package edu.mit.broad.picard.sam;

import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.metrics.MetricsFile;
import edu.mit.broad.picard.util.Log;
import edu.mit.broad.sam.util.SortingCollection;
import edu.mit.broad.sam.*;

import java.io.*;
import java.util.*;

/**
 * A better duplication marking algorithm that handles all cases including clipped
 * and gapped alignments.
 *
 * @author Tim Fennell
 */
public class MarkDuplicates extends CommandLineProgram {
    private static final Log log = Log.getInstance(MarkDuplicates.class);

    @Usage public final String USAGE =
            "Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules. " +
            "All records are then written to the output file with the duplicate records flagged.";
    @Option(shortName="I", doc="The input SAM or BAM file to analyze") public File INPUT;
    @Option(shortName="O", doc="The output file to right marked records to") public File OUTPUT;
    @Option(shortName="M", doc="File to write duplication metrics to") public File METRICS_FILE;

    private SortingCollection<ReadEnds> pairSort;
    private SortingCollection<ReadEnds> fragSort;
    private long[] duplicateIndexes = new long[1000000];
    private int nextIndex = 0; // The next offset into duplicateIndexes to use


    /** Stock main method. */
    public static void main(String[] args) {
        new MarkDuplicates().instanceMain(args);
    }

    /** Little struct-like class to hold read pair (and fragment) end data. */
    private static class ReadEnds {
        public static final int SIZE_OF = (1*1) + (2*1) + (4*4) + (8*2) + 8; // last 8 == reference overhead
        public static final byte F=0, R=1, FF=2, FR=3, RR=4, RF=5;

        short score = 0;
        byte orientation;
        int read1Sequence     = -1;
        int read1Coordinate   = -1;
        long read1IndexInFile = -1;
        int read2Sequence     = -1;
        int read2Coordinate   = -1;
        long read2IndexInFile = -1;

        boolean isPaired() { return this.read2Sequence != -1; }
    }

    /** Comparator for ReadEnds that orders by read1 position then pair orientation then read2 position. */
    private static class ReadEndsComparator implements Comparator<ReadEnds> {
        public int compare(ReadEnds lhs, ReadEnds rhs) {
            int retval = lhs.read1Sequence - rhs.read1Sequence;
            if (retval == 0) retval = lhs.read1Coordinate - rhs.read1Coordinate;
            if (retval == 0) retval = lhs.orientation - rhs.orientation;
            if (retval == 0) retval = lhs.read2Sequence   - rhs.read2Sequence;
            if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
            if (retval == 0) retval = (int) (lhs.read1IndexInFile - rhs.read1IndexInFile);
            if (retval == 0) retval = (int) (lhs.read2IndexInFile - rhs.read2IndexInFile);

            return retval;
        }
    }

    /** Coded for ReadEnds that just outputs the primitive fields and reads them back. */
    private static class ReadEndsCodec implements SortingCollection.Codec<ReadEnds> {
        private DataInputStream in;
        private DataOutputStream out;

        public SortingCollection.Codec<ReadEnds> clone() {
            return new ReadEndsCodec();
        }

        public void setOutputStream(OutputStream os) { this.out = new DataOutputStream(os); }
        public void setInputStream(InputStream is) { this.in = new DataInputStream(is); }

        public void encode(ReadEnds read) {
            try {
                this.out.writeShort(read.score);
                this.out.writeByte(read.orientation);
                this.out.writeInt(read.read1Sequence);
                this.out.writeInt(read.read1Coordinate);
                this.out.writeLong(read.read1IndexInFile);
                this.out.writeInt(read.read2Sequence);

                if (read.orientation > ReadEnds.R) {
                    this.out.writeInt(read.read2Coordinate);
                    this.out.writeLong(read.read2IndexInFile);
                }
                this.out.flush();
            }
            catch (IOException ioe) {
                throw new PicardException("Exception writing ReadEnds to file.", ioe);
            }
        }

        public ReadEnds decode() {
            ReadEnds read = new ReadEnds();
            try {
                // If the first read results in an EOF we've exhausted the stream
                try { read.score = this.in.readShort(); }
                catch (EOFException eof) { return null; }

                read.orientation      = this.in.readByte();
                read.read1Sequence    = this.in.readInt();
                read.read1Coordinate  = this.in.readInt();
                read.read1IndexInFile = this.in.readLong();
                read.read2Sequence    = this.in.readInt();

                if (read.orientation > ReadEnds.R) {
                    read.read2Coordinate  = this.in.readInt();
                    read.read2IndexInFile = this.in.readLong();
                }
                return read;
            }
            catch (IOException ioe) {
                throw new PicardException("Exception writing ReadEnds to file.", ioe);
            }
        }
    }

    /**
     * Main work method.  Reads the BAM file once and collects sorted information about
     * the 5' ends of both ends of each read (or just one end in the case of pairs).
     * Then makes a pass through those determining duplicates before re-reading the
     * input file and writing it out with duplication flags set correctly.
     */
    protected int doWork() {
        log.info("Reading input file and constructing read end information.");
        buildSortedReadEndLists();
        generateDuplicateIndexes();
        log.info("Marking " + this.duplicateIndexes.length + " records as duplicates.");
        DuplicationMetrics metrics = new DuplicationMetrics();
        SAMFileReader in  = new SAMFileReader(INPUT);
        SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(),
                                                                          true,
                                                                          OUTPUT);

        // Now copy over the file while marking all the necessary indexes as duplicates
        long recordInFileIndex = 0;
        long nextDuplicateIndex = (this.duplicateIndexes.length == 0 ? -1 :  this.duplicateIndexes[0]);
        int  arrayIndex = 1;

        for (SAMRecord rec : in) {
            // First bring the simple metrics up to date
            if (rec.getReadUnmappedFlag()) {
                ++metrics.UNMAPPED_READS;
            }
            else if (!rec.getReadPairedFlag() || rec.getMateUnmappedFlag()) {
                ++metrics.UNPAIRED_READS_EXAMINED;
            }
            else if (rec.getFirstOfPairFlag()){
                ++metrics.READ_PAIRS_EXAMINED;
            }


            if (recordInFileIndex++ == nextDuplicateIndex) {
                rec.setDuplicateReadFlag(true);

                // Update the duplication metrics
                if (!rec.getReadPairedFlag() || rec.getMateUnmappedFlag()) {
                    ++metrics.UNPAIRED_READ_DUPLICATES;
                }
                else if (rec.getFirstOfPairFlag()) {
                    ++metrics.READ_PAIR_DUPLICATES;
                }

                // Now try and figure out the next duplicate index
                try {
                    nextDuplicateIndex = this.duplicateIndexes[arrayIndex++];
                }
                catch (ArrayIndexOutOfBoundsException e) {
                    // Only happens once we've marked all the duplicates
                    nextDuplicateIndex = -1;
                    arrayIndex = -1;
                }
            }

            out.addAlignment(rec);
        }

        out.close();


        // Write out the metrics
        metrics.calculateDerivedMetrics();
        MetricsFile<DuplicationMetrics,Double> file = getMetricsFile();
        file.addMetric(metrics);
        file.setHistogram(metrics.calculateRoiHistogram());
        file.write(METRICS_FILE);        

        return 0;
    }

    /**
     * Goes through all the records in a file and generates a set of ReadEnds objects that
     * hold the necessary information (reference sequence, 5' read coordinate) to do
     * duplication, caching to disk as necssary to sort them.
     */
    private void buildSortedReadEndLists() {
        // TODO: take into account clipping/padding?
        int maxInMemory = (int) ((Runtime.getRuntime().maxMemory() * 0.25) / ReadEnds.SIZE_OF);
        this.pairSort = SortingCollection.newInstance(ReadEnds.class,
                                                      new ReadEndsCodec(),
                                                      new ReadEndsComparator(),
                                                      maxInMemory);

        this.fragSort = SortingCollection.newInstance(ReadEnds.class,
                                                      new ReadEndsCodec(),
                                                      new ReadEndsComparator(),
                                                      maxInMemory);

        Map<String, ReadEnds> tmp = new HashMap<String, ReadEnds>();
        SAMFileReader sam = new SAMFileReader(INPUT);
        SAMFileHeader header = sam.getFileHeader();
        long index = 0;

        for (SAMRecord rec : sam) {
            if (rec.getReadUnmappedFlag()) {
                continue;
            }

            ReadEnds fragmentEnd = buildReadEnds(header, index, rec);
            this.fragSort.add(fragmentEnd);

            if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
                String key = rec.getAttribute(ReservedTagConstants.READ_GROUP_ID) + ":" + rec.getReadName();
                ReadEnds pairedEnds = tmp.remove(key);

                // See if we've already seen the first end or not
                if (pairedEnds == null) {
                    pairedEnds = buildReadEnds(header, index, rec);
                    tmp.put(key, pairedEnds);
                }
                else {
                    int sequence = fragmentEnd.read1Sequence;
                    int coordinate = fragmentEnd.read1Coordinate;

                    // If the second read is actually later, just add the second read data, else flip the reads
                    if (sequence > pairedEnds.read1Sequence || (sequence == pairedEnds.read1Sequence && coordinate >= pairedEnds.read1Coordinate)) {
                        pairedEnds.read2Sequence    = sequence;
                        pairedEnds.read2Coordinate  = coordinate;
                        pairedEnds.read2IndexInFile = index;
                        pairedEnds.orientation = getOrientationByte(pairedEnds.orientation == ReadEnds.R, rec.getReadNegativeStrandFlag());
                    }
                    else {
                        pairedEnds.read2Sequence    = pairedEnds.read1Sequence;
                        pairedEnds.read2Coordinate  = pairedEnds.read1Coordinate;
                        pairedEnds.read2IndexInFile = pairedEnds.read1IndexInFile;
                        pairedEnds.read1Sequence    = sequence;
                        pairedEnds.read1Coordinate  = coordinate;
                        pairedEnds.read1IndexInFile = index;
                        pairedEnds.orientation = getOrientationByte(rec.getReadNegativeStrandFlag(), pairedEnds.orientation == ReadEnds.R);
                    }

                    pairedEnds.score += getScore(rec);
                    this.pairSort.add(pairedEnds);
                }
            }

            ++index;
        }
    }

    /** Builds a read ends object that represents a single read. */
    private ReadEnds buildReadEnds(SAMFileHeader header, long index, SAMRecord rec) {
        ReadEnds ends = new ReadEnds();
        ends.read1Sequence    = rec.getReferenceIndex(header);
        ends.read1Coordinate  = rec.getReadNegativeStrandFlag() ? rec.getUnclippedEnd() : rec.getUnclippedStart();
        ends.orientation = rec.getReadNegativeStrandFlag() ? ReadEnds.R : ReadEnds.F;
        ends.read1IndexInFile = index;
        ends.score = getScore(rec);

        // Doing this lets the ends object know that it's part of a pair
        if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
            ends.read2Sequence = rec.getMateReferenceIndex(header);
        }

        return ends;
    }

    /**
     * Returns a single byte that encodes the orientation of the two reads in a pair.
     */
    private byte getOrientationByte(boolean read1NegativeStrand, boolean read2NegativeStrand) {
        if (read1NegativeStrand) {
            if (read2NegativeStrand)  return ReadEnds.RR;
            else return ReadEnds.RF;
        }
        else {
            if (read2NegativeStrand)  return ReadEnds.FR;
            else return ReadEnds.FF;
        }
    }



    /** Calculates a score for the read which is the sum of scores over Q20. */
    private short getScore(SAMRecord rec) {
        short score = 0;
        for (byte b : rec.getBaseQualities()) {
            if (b >= 15) score += b;
        }

        return score;
    }

    /**
     * Goes through the accumulated ReadEnds objects and determines which of them are
     * to be marked as duplicates.
     *
     * @return an array with an ordered list of indexes into the source file
     */
    private void generateDuplicateIndexes() {
        ReadEnds firstOfNextChunk = null;
        List<ReadEnds> nextChunk  = new ArrayList<ReadEnds>(200);

        // First just do the pairs
        log.info("Traversing read pair information and detecting duplicates.");
        for (ReadEnds next : this.pairSort) {
            if (firstOfNextChunk == null) {
                firstOfNextChunk = next;
                nextChunk.add(firstOfNextChunk);
            }
            else if (areComparableForDuplicates(firstOfNextChunk, next, true)) {
                nextChunk.add(next);
            }
            else {
                if (nextChunk.size() > 1) {
                    markDuplicatePairs(nextChunk);
                }

                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
            }
        }
        markDuplicatePairs(nextChunk);
        this.pairSort = null;

        // Now deal with the fragments
        log.info("Traversing fragment information and detecting duplicates.");
        boolean containsPairs = false;
        boolean containsFrags = false;

        for (ReadEnds next : this.fragSort) {
            if (firstOfNextChunk != null && areComparableForDuplicates(firstOfNextChunk, next, false)) {
                nextChunk.add(next);
                containsPairs = containsPairs || next.isPaired();
                containsFrags = containsFrags || !next.isPaired();
            }
            else {
                if (nextChunk.size() > 1 && containsFrags) {
                    markDuplicateFragments(nextChunk, containsPairs);
                }

                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
                containsPairs = next.isPaired();
                containsFrags = !next.isPaired();
            }
        }
        markDuplicateFragments(nextChunk, containsPairs);
        this.fragSort = null;

        // Now shrink down the array and sort it
        log.info("Sorting list of duplicate records.");
        long[] tmp = new long[this.nextIndex];
        System.arraycopy(this.duplicateIndexes, 0, tmp, 0, this.nextIndex);
        this.duplicateIndexes = tmp;
        Arrays.sort(this.duplicateIndexes);
    }

    private boolean areComparableForDuplicates(final ReadEnds lhs, final ReadEnds rhs, final boolean compareRead2) {
        boolean retval =  lhs.read1Sequence   == rhs.read1Sequence &&
                          lhs.read1Coordinate == rhs.read1Coordinate &&
                          lhs.orientation     == rhs.orientation;

        if (compareRead2) {
            retval = lhs.read2Sequence   == rhs.read2Sequence &&
                     lhs.read2Coordinate == rhs.read2Coordinate;
        }

        return retval;
    }

    private void addIndexAsDuplicate(final long bamIndex) {
        if (this.nextIndex > this.duplicateIndexes.length - 1) {
            long[] tmp = new long[this.duplicateIndexes.length * 2];
            System.arraycopy(this.duplicateIndexes, 0, tmp, 0, this.nextIndex);
            this.duplicateIndexes = tmp;
        }

        this.duplicateIndexes[this.nextIndex++] = bamIndex;
    }

    /**
     * Takes a list of ReadEnds objects and removes from it all objects that should
     * not be marked as duplicates.
     *
     * @param list
     */
    private void markDuplicatePairs(final List<ReadEnds> list) {
        short maxScore = 0;
        ReadEnds best = null;

        for (final ReadEnds end : list) {
            if (end.score > maxScore || best == null) {
                maxScore = end.score;
                best = end;
            }
        }

        for (final ReadEnds end : list) {
            if (end != best) {
                addIndexAsDuplicate(end.read1IndexInFile);
                addIndexAsDuplicate(end.read2IndexInFile);
            }
        }
    }

    /**
     * Takes a list of ReadEnds objects and removes from it all objects that should
     * not be marked as duplicates.
     *
     * @param list
     */
    private void markDuplicateFragments(final List<ReadEnds> list, final boolean containsPairs) {
        if (containsPairs) {
            for (final ReadEnds end : list) {
                if (!end.isPaired()) addIndexAsDuplicate(end.read1IndexInFile);
            }
        }
        else {
            short maxScore = 0;
            ReadEnds best = null;
            for (final ReadEnds end : list) {
                if (end.score > maxScore || best == null) {
                    maxScore = end.score;
                    best = end;
                }
            }

            for (final ReadEnds end : list) {
                if (end != best) {
                    addIndexAsDuplicate(end.read1IndexInFile);
                }
            }
        }
    }
}
