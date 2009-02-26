package edu.mit.broad.picard.quality;

import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.sam.SAMRecord;
import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.sam.AlignmentBlock;
import edu.mit.broad.sam.util.CloseableIterator;
import edu.mit.broad.picard.variation.DbSnpFileReader;
import edu.mit.broad.picard.variation.KnownVariant;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequence;
import edu.mit.broad.picard.util.CoordMath;
import edu.mit.broad.picard.util.Histogram;
import edu.mit.broad.picard.util.SequenceUtil;

import java.util.Map;
import java.util.BitSet;
import java.util.TreeMap;

/**
 * Takes a set of aligned reads with qualities and determines the empirical quality
 * score for each of the bins.
 *
 * @author Tim Fennell
 */
public class QualityScoreCalibrator {
    private final SAMFileReader sam;
    private final ReferenceSequenceFile ref;
    private final DbSnpFileReader dbsnp;

    private QualityScoreMatrix read1Matrix;
    private QualityScoreMatrix read2Matrix;

    /**
     * Constructs a calibrator that will read records from the specified SAMFileReader
     * and compare them the supplied reference. Optionally takes a set of known variants
     * who's positions will be excluded during calibration.
     *
     * @param sam the set of SAM records to use to calibrate qualities
     * @param ref the reference sequence against which the records were aligned
     * @param dbsnp the (optional) set of dbsnp positions to mask during calibration
     */
    public QualityScoreCalibrator(SAMFileReader sam, ReferenceSequenceFile ref, DbSnpFileReader dbsnp) {
        this.sam = sam;
        this.dbsnp = dbsnp;
        this.ref = ref;
    }

    /**
     * Calculates calibrated quality scores using at most the specified number of aligned
     * reads. If the end of the file is hit first then fewer reads will be used.
     *
     * @param readLimit the number of aligned reads to use if the file contains more
     */
    public void calibrate(final int readLimit) {
        ReferenceSequence reference = null;
        SAMFileHeader header = this.sam.getFileHeader();
        CloseableIterator<SAMRecord> samIterator = this.sam.iterator();
        SAMRecord read = samIterator.next();
        int readsProcessed = 0;

        // Quality score matrixes for reads 1 and 2 separately
        this.read1Matrix = new QualityScoreMatrix();
        this.read2Matrix = new QualityScoreMatrix();


        refloop: while ((reference = this.ref.nextSequence()) != null) {
            final byte[] refBases = reference.getBases();
            final BitSet snps = getDbSnpMask(reference);

            while (read != null && read.getReferenceIndex(header) == reference.getContigIndex()) {
                if (!read.getReadUnmappedFlag() && !read.getNotPrimaryAlignmentFlag()) {
                    final QualityScoreMatrix matrix = read.getFirstOfPairFlag() ? this.read1Matrix : this.read2Matrix;
                    final byte[] readBases = read.getReadBases();
                    final byte[] qualities  = read.getBaseQualities();

                    for (AlignmentBlock block : read.getAlignmentBlocks()) {
                        final int readIndex = block.getReadStart() - 1;
                        final int refIndex  = block.getReferenceStart() - 1;
                        final int length    = block.getLength();

                        for (int i=0; i<length; ++i) {
                            // Skip dbSNP loci
                            if (snps.get(refIndex+i+1)) continue;

                            final int readBaseIndex = readIndex+i;
                            boolean match = SequenceUtil.basesEqual(readBases[readBaseIndex], refBases[refIndex+i]);
                            int cycle = CoordMath.getCycle(
                                    read.getReadNegativeStrandFlag(), readBases.length, readBaseIndex); 
                            matrix.addObservation(cycle, qualities[readBaseIndex], !match);
                        }
                    }

                    if (readLimit > 0 && ++readsProcessed >= readLimit) {
                        break refloop;
                    }
                }

                // Advance the sam iterator
                if (samIterator.hasNext()) {
                    read = samIterator.next();
                }
                else {
                    read = null;
                }
            }
        }

        this.read1Matrix.computeCalibratedQualities();
        if (!this.read2Matrix.isEmpty()) this.read2Matrix.computeCalibratedQualities();
    }

    /** Gets the calibration matrix for the first read. */
    public QualityScoreMatrix getRead1Matrix() { return read1Matrix; }

    /** Gets the calibration matrix for the second read. May be empty if there was no second read data. */
    public QualityScoreMatrix getRead2Matrix() { return read2Matrix; }

    /**
     * Returns a BitSet that denotes whether a dbSNP entry is present at each
     * base in the reference sequence.  The set is reference.length() + 1 so that
     * it can be indexed by 1-based reference base.  True means dbSNP present,
     * false means no dbSNP present.
     */
    private BitSet getDbSnpMask(ReferenceSequence reference) {
        int index = reference.getContigIndex();
        BitSet bits = new BitSet(reference.length() + 1);

        /* Just return an all false bit set if we don't have dbsnp data. */
        if (this.dbsnp == null) {
            return bits;
        }

        /* Read off the next contig's worth of data. */
        while (this.dbsnp.hasNext()) {
            KnownVariant variant = this.dbsnp.peek();

            if (variant.getSequenceIndex() < index) {
                this.dbsnp.next();
            }
            else if (variant.getSequenceIndex() == index) {
                variant = this.dbsnp.next();

                for (int i=variant.getStartPos(); i<=variant.getEndPos(); ++i) {
                    bits.set(i, true);
                }
            }
            else {
                break;
            }
        }

        return bits;
    }
}
