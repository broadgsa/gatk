/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.picard.sam;

import java.io.File;

import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.picard.metrics.AggregateMetricCollector;
import edu.mit.broad.picard.metrics.MetricBase;
import edu.mit.broad.picard.metrics.MetricCollector;
import edu.mit.broad.picard.metrics.MetricsFile;
import edu.mit.broad.picard.metrics.StringHeader;
import edu.mit.broad.picard.reference.ReferenceSequence;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;
import edu.mit.broad.picard.sam.CollectAlignmentSummaryMetrics.AlignmentSummaryMetrics.Type;
import edu.mit.broad.picard.util.CoordMath;
import edu.mit.broad.picard.util.Histogram;
import edu.mit.broad.picard.util.SequenceUtil;
import edu.mit.broad.sam.AlignmentBlock;
import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.sam.SAMRecord;
import edu.mit.broad.sam.util.CloseableIterator;

/**
 * A command line tool to read a BAM file and produce standard alignment metrics that would be applicable to any alignment.  
 * Metrics to include, but not limited to:
 * <ul>
 * <li>Total number of reads (total, period, no exclusions)</li>
 * <li>Total number of PF reads (PF == does not fail vendor check flag)</li>
 * <li>Number of PF noise reads (does not fail vendor check and has noise attr set)</li>
 * <li>Total aligned PF reads (any PF read that has a sequence and position)</li>
 * <li>High quality aligned PF reads (high quality == mapping quality >= 20)</li>
 * <li>High quality aligned PF bases (actual aligned bases, calculate off alignment blocks)</li>
 * <li>High quality aligned PF Q20 bases (subset of above where base quality >= 20)</li>
 * <li>Median mismatches in HQ aligned PF reads (how many aligned bases != ref on average)</li>
 * <li>Reads aligned in pairs (vs. reads aligned with mate unaligned/not present)</li>
 * <li>Read length (how to handle mixed lengths?)</li>
 * <li>Bad Cycles - how many machine cycles yielded combined no-call and mismatch rates of >= 80%</li>
 * <li>Strand balance - reads mapped to positive strand / total mapped reads</li>
 * </ul>
 * Metrics are written for the first read of a pair, the second read, and combined for the pair.
 * 
 * @author Doug Voet
 */
public class CollectAlignmentSummaryMetrics extends CommandLineProgram {
    private static final int MAPPING_QUALITY_THRESHOLD = 20;
    private static final int BASE_QUALITY_THRESHOLD = 20;

    // Usage and parameters
    @Usage(programVersion="1.0") 
    public String USAGE = "Reads a SAM or BAM file and writes a file containing summary metrics.\n";
    @Option(shortName="I", doc="SAM or BAM file") public File INPUT;
    @Option(shortName="O", doc="File to write insert size metrics to") public File OUTPUT;
    @Option(shortName="R", doc="Reference sequence file") public File REFERENCE;
    @Option(doc="If true (default), \"unsorted\" SAM/BAM files will be considerd coordinate sorted")
    public Boolean ASSUME_COODINATE_SORTED = Boolean.TRUE;

    private ReferenceSequenceFile ref;
    private ReferenceSequence refSequence;
    private SAMFileHeader samFileHeader;

    /** Required main method implementation. */
    public static void main(String[] argv) {
        System.exit(new CollectAlignmentSummaryMetrics().instanceMain(argv));
    }

    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsReadable(REFERENCE);
        IoUtil.assertFileIsWritable(OUTPUT);
        SAMFileReader in = new SAMFileReader(INPUT);
        assertCoordinateSortOrder(in);

        this.ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE);
        this.samFileHeader = in.getFileHeader();

        MetricsFile<AlignmentSummaryMetrics, Comparable<?>> file = collectMetrics(in.iterator());
        in.close();

        file.write(OUTPUT);

        return 0;
    }

    private void assertCoordinateSortOrder(SAMFileReader in) {
        switch (in.getFileHeader().getSortOrder()) {
        case coordinate:
            break;
        case unsorted:
            if (this.ASSUME_COODINATE_SORTED) {
                break;
            }
        default:
            throw new PicardException("Cannot collect summary statistics in file " + INPUT.getAbsoluteFile() +
            " because it is not sorted in coordinate order.");
        }
    }

    private ReferenceSequence getReference(SAMRecord record) {
        while (refSequence == null || 
                record.getReferenceIndex(samFileHeader) > refSequence.getContigIndex()) {

            refSequence = ref.nextSequence();
        }

        if (refSequence == null || record.getReferenceIndex() != refSequence.getContigIndex()) {
            throw new PicardException("Cannot find reference sequence [" + 
                    record.getReferenceIndex() + "] in reference file");
        }

        return refSequence;
    }

    /**
     * Does all the work of iterating through the sam file and collecting summary alignment metrics.
     */
    private MetricsFile<AlignmentSummaryMetrics, Comparable<?>> collectMetrics(
            CloseableIterator<SAMRecord> samIterator) {

        final MetricCollector<AlignmentSummaryMetrics> unpairedCollector = 
            constructCollector(Type.UNPAIRED);
        final MetricCollector<AlignmentSummaryMetrics> firstOfPairCollector = 
            constructCollector(Type.FIRST_OF_PAIR);
        final MetricCollector<AlignmentSummaryMetrics> secondOfPairCollector = 
            constructCollector(Type.SECOND_OF_PAIR);
        final MetricCollector<AlignmentSummaryMetrics> pairCollector = 
            constructCollector(Type.PAIR);

        while (samIterator.hasNext()) {
            SAMRecord record = samIterator.next();

            if (record.getReadPairedFlag()) {
                if (record.getFirstOfPairFlag()) {
                    firstOfPairCollector.addRecord(record);
                } else {
                    secondOfPairCollector.addRecord(record);
                }
                pairCollector.addRecord(record);
            } else {
                unpairedCollector.addRecord(record);
            }
        }

        firstOfPairCollector.onComplete();
        secondOfPairCollector.onComplete();
        pairCollector.onComplete();
        unpairedCollector.onComplete();
        
        MetricsFile<AlignmentSummaryMetrics, Comparable<?>> file = getMetricsFile();
        file.addHeader(new StringHeader("Input file: " + INPUT.getAbsolutePath()));
        file.addHeader(new StringHeader("Output file: " + OUTPUT.getAbsolutePath()));
        file.addHeader(new StringHeader("Reference file: " + REFERENCE.getAbsolutePath()));
        
        if (firstOfPairCollector.getMetrics().TOTAL_READS > 0) {
            file.addMetric(firstOfPairCollector.getMetrics());
            // override how bad cycle is determined for paired reads, it should be
            // the sum of first and second reads
            pairCollector.getMetrics().BAD_CYCLES = 
                firstOfPairCollector.getMetrics().BAD_CYCLES +
                secondOfPairCollector.getMetrics().BAD_CYCLES;
            file.addMetric(secondOfPairCollector.getMetrics());
            file.addMetric(pairCollector.getMetrics());
        }
        if (unpairedCollector.getMetrics().TOTAL_READS > 0) {
            file.addMetric(unpairedCollector.getMetrics());
        }

        return file;
    }

    private MetricCollector<AlignmentSummaryMetrics> constructCollector(Type type) {
        MetricCollector<AlignmentSummaryMetrics> collector = 
            new AggregateMetricCollector<AlignmentSummaryMetrics>(new ReadCounter(), new QualityMappingCounter());
        collector.setMetrics(new AlignmentSummaryMetrics());
        collector.getMetrics().TYPE = type;
        return collector;
    }

    public static class AlignmentSummaryMetrics extends MetricBase {
        public enum Type { UNPAIRED, FIRST_OF_PAIR, SECOND_OF_PAIR, PAIR }
        public Type TYPE;
        public long TOTAL_READS;
        public long PF_READS;
        public long PF_NOISE_READS;
        public long PF_READS_ALIGNED;
        public long PF_HQ_ALIGNED_READS;
        public long PF_HQ_ALIGNED_BASES;
        public long PF_HQ_ALIGNED_Q20_BASES;
        public double PF_HQ_MEDIAN_MISMATCHES;
        public double MEAN_READ_LENGTH;
        public long READS_ALIGNED_IN_PAIRS;
        public long BAD_CYCLES;
        public double STRAND_BALANCE;
    }

    /** counts reads that match various conditions */
    private class ReadCounter implements MetricCollector<AlignmentSummaryMetrics> {
        private long numPositiveStrand = 0;
        private Histogram<Integer> readLengthHistogram = new Histogram<Integer>();
        private AlignmentSummaryMetrics metrics;

        @Override
        public void addRecord(SAMRecord record) {
            if (record.getNotPrimaryAlignmentFlag()) {
                // only want 1 count per read so skip non primary alignments
                return;
            }

            metrics.TOTAL_READS++;
            readLengthHistogram.increment(record.getReadBases().length);

            if (!record.getReadFailsVendorQualityCheckFlag()) {
                metrics.PF_READS++;

                if (isNoiseRead(record)) {
                    metrics.PF_NOISE_READS++;
                }
                if (!record.getReadUnmappedFlag()) {
                    metrics.PF_READS_ALIGNED++;
                }
            }

            if (!record.getReadUnmappedFlag() && 
                    record.getReadPairedFlag() &&
                    !record.getMateUnmappedFlag()) {
                metrics.READS_ALIGNED_IN_PAIRS++;
            }

            if (!record.getReadNegativeStrandFlag()) {
                numPositiveStrand++;
            }
        }

        @Override
        public void onComplete() {
            metrics.MEAN_READ_LENGTH = readLengthHistogram.getMean();
            metrics.STRAND_BALANCE = numPositiveStrand / (double) metrics.TOTAL_READS;
        }

        private boolean isNoiseRead(SAMRecord record) {
            final Object noiseAttribute = record.getAttribute(ReservedTagConstants.XN);
            return (noiseAttribute != null && noiseAttribute.equals(1));
        }

        @Override
        public void setMetrics(AlignmentSummaryMetrics metrics) {
            this.metrics = metrics;
        }

        @Override
        public AlignmentSummaryMetrics getMetrics() {
            return this.metrics;
        }
    }

    /** counts quality mappings & base calls that match various conditions */
    private class QualityMappingCounter implements MetricCollector<AlignmentSummaryMetrics> {
        private Histogram<Long> mismatchHistogram = new Histogram<Long>();
        private Histogram<Integer> badCycleHistogram = new Histogram<Integer>();
        private AlignmentSummaryMetrics metrics;

        @Override
        public void addRecord(SAMRecord record) {
            if (record.getNotPrimaryAlignmentFlag()) {
                return;
            }
            if (record.getReadUnmappedFlag()) {
                final byte[] readBases = record.getReadBases();
                for (int i = 0; i < readBases.length; i++) {
                    if (SequenceUtil.isNoCall(readBases[i])) {
                        badCycleHistogram.increment(CoordMath.getCycle(record.getReadNegativeStrandFlag(), readBases.length, i));
                    }
                }
            } else {
                boolean highQualityMapping = isHighQualityMapping(record);
                if (highQualityMapping) metrics.PF_HQ_ALIGNED_READS++;
                
                final byte[] readBases = record.getReadBases();
                final byte[] refBases = getReference(record).getBases();
                final byte[] qualities  = record.getBaseQualities();
                long mismatchCount = 0;
                
                for (AlignmentBlock alignmentBlock : record.getAlignmentBlocks()) {
                    final int readIndex = alignmentBlock.getReadStart() - 1;
                    final int refIndex  = alignmentBlock.getReferenceStart() - 1;
                    final int length    = alignmentBlock.getLength();
                    if (highQualityMapping) metrics.PF_HQ_ALIGNED_BASES += alignmentBlock.getLength();
    
                    for (int i=0; i<length; ++i) {
                        final int readBaseIndex = readIndex + i;
                        boolean mismatch = !SequenceUtil.basesEqual(readBases[readBaseIndex], refBases[refIndex+i]);
                        if (highQualityMapping) {
                            if (qualities[readBaseIndex] >= BASE_QUALITY_THRESHOLD) {
                                metrics.PF_HQ_ALIGNED_Q20_BASES++;
                            }
                            if (mismatch) {
                                mismatchCount++;
                            }
                        }
                        if (mismatch || SequenceUtil.isNoCall(readBases[readBaseIndex])) {
                            badCycleHistogram.increment(CoordMath.getCycle(record.getReadNegativeStrandFlag(), readBases.length, i));
                        }
                    }
                }
                mismatchHistogram.increment(mismatchCount);
            }
        }

        private boolean isHighQualityMapping(SAMRecord record) {
            return !record.getReadFailsVendorQualityCheckFlag() &&
            record.getMappingQuality() >= MAPPING_QUALITY_THRESHOLD;
        }

        @Override
        public void onComplete() {
            metrics.PF_HQ_MEDIAN_MISMATCHES = mismatchHistogram.getMedian();
            metrics.BAD_CYCLES = 0;

            for (Histogram<Integer>.Bin cycleBin : badCycleHistogram.values()) {
                double badCyclePercentage = cycleBin.getValue() / metrics.TOTAL_READS;
                if (badCyclePercentage >= .8) {
                    metrics.BAD_CYCLES++;
                }
            }
        }

        @Override
        public void setMetrics(AlignmentSummaryMetrics metrics) {
            this.metrics = metrics;
        }

        @Override
        public AlignmentSummaryMetrics getMetrics() {
            return this.metrics;
        }
    }
}
