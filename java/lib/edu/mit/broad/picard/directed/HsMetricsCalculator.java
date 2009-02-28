package edu.mit.broad.picard.directed;

import edu.mit.broad.picard.util.*;
import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.sam.SAMRecord;
import edu.mit.broad.sam.AlignmentBlock;
import edu.mit.broad.sam.SAMSequenceRecord;

import java.util.*;
import java.io.*;

/**
 * Calculates HS metrics for a given SAM or BAM file. Requires the input of a list of
 * target intervals and a list of bait intervals. Can be invoked either on an entire
 * iterator of SAMRecords or be passed SAMRecords one at a time.
 *
 * @author Tim Fennell
 */
public class HsMetricsCalculator {
    // What is considered "near" to the bait
    private static final int NEAR_BAIT_DISTANCE = 250;
    private static final Log log = Log.getInstance(HsMetricsCalculator.class);

    // Holds file names and other parameter related junk
    private SAMFileReader sam;
    private File baitFile;
    private File targetFile;
    private IntervalList baits;
    private IntervalList targets;

    // Overlap detector for finding overlaps between reads and the experimental targets
    private OverlapDetector<Interval> targetDetector = new OverlapDetector<Interval>(0,0);

	// Overlap detector for finding overlaps between the reads and the baits (and the near bait space)
    private OverlapDetector<Interval> baitDetector = new OverlapDetector<Interval>(-NEAR_BAIT_DISTANCE,0);

    // A Map to accumulate per-bait-region (i.e. merge of overlapping baits) coverage. */
    private Map<Interval, Coverage> coverageByTarget = new HashMap<Interval, Coverage>();

    private HsMetrics metrics = new HsMetrics();

    /**
	 * Constructor that parses the squashed reference to genome reference file and stores the
	 * information in a map for later use.
	 */
    public HsMetricsCalculator(File baits, File targets) {
        this.baitFile     = baits;
        this.targetFile   = targets;
        this.baits        = IntervalList.fromFile(baits);
        this.targets      = IntervalList.fromFile(targets);

        this.metrics.BAIT_SET = baits.getName();
        int tmp = this.metrics.BAIT_SET.indexOf(".");
        if (tmp > 0) {
            this.metrics.BAIT_SET = this.metrics.BAIT_SET.substring(0, tmp);
        }

        List<Interval> uniqueBaits = this.baits.getUniqueIntervals();
        this.baitDetector.addAll(uniqueBaits, uniqueBaits);
        this.metrics.BAIT_TERRITORY = Interval.countBases(uniqueBaits);

        List<Interval> uniqueTargets = this.targets.getUniqueIntervals();
        this.targetDetector.addAll(uniqueTargets, uniqueTargets);
        this.metrics.TARGET_TERRITORY = Interval.countBases(uniqueTargets);

        for (SAMSequenceRecord seq : this.baits.getHeader().getSequences()) {
            this.metrics.GENOME_SIZE += seq.getSequenceLength();
        }

        // Populate the coverage by target map
        for (Interval target : this.targets.getIntervals()) {
            this.coverageByTarget.put(target, new Coverage(target, 0));
        }
    }

    /** Iterates over all records in the file and collects metrics. */
    public void analyze(Iterator<SAMRecord> records) {
        int i = 0;
        while (records.hasNext()) {
            analyze(records.next());

            if (++i % 1000000 == 0) {
                log.info("Processed " + i + " records so far.");
            }
        }
    }

    /** Adds information about an individual SAMRecord to the statistics. */
    public void analyze(SAMRecord rec) {
        // Just plain avoid records that are marked as not-primary
        if (rec.getNotPrimaryAlignmentFlag()) return;
        
        this.metrics.TOTAL_READS += 1;

        // Check for PF reads
        if (rec.getReadFailsVendorQualityCheckFlag()) {
            return;
        }
        else {
            ++this.metrics.PF_READS;
        }

        // Check for reads that are marked as duplicates
        if (rec.getDuplicateReadFlag()) {
            return;
        }
        else {
            ++this.metrics.PF_UNIQUE_READS;
        }

        // Don't bother with reads that didn't align uniquely
        if (rec.getReadUnmappedFlag() || rec.getMappingQuality() == 0) {
            return;
        }

        this.metrics.PF_READS_ALIGNED += 1;
        for (AlignmentBlock block : rec.getAlignmentBlocks()) {
            this.metrics.PF_BASES_ALIGNED += block.getLength();
        }

        Interval read = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());

        // Find the target overlaps
        Collection<Interval> targets = this.targetDetector.getOverlaps(read);
        if (targets != null && !targets.isEmpty()) {
            for (Interval target : targets) {
                Coverage coverage = this.coverageByTarget.get(target);

                for (AlignmentBlock block : rec.getAlignmentBlocks()) {
                    int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());
                    for (int pos=block.getReferenceStart(); pos<=end; ++ pos) {
                        if (pos >= target.getStart() && pos <= target.getEnd()) {
                            ++this.metrics.ON_TARGET_BASES;
                            coverage.addBase(pos - target.getStart());
                        }
                    }
                }
            }
        }

        // Now do the bait overlaps
        int mappedBases = 0;
        for (AlignmentBlock block : rec.getAlignmentBlocks()) mappedBases += block.getLength();
        Collection<Interval> baits = this.baitDetector.getOverlaps(read);
        int onBaitBases = 0;

        if (baits != null && !baits.isEmpty()) {
            for (Interval bait : baits) {
                for (AlignmentBlock block : rec.getAlignmentBlocks()) {
                    int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());

                    for (int pos=block.getReferenceStart(); pos<=end; ++pos) {
                        if (pos >= bait.getStart() && pos <= bait.getEnd()) ++onBaitBases;
                    }
                }
            }

            this.metrics.ON_BAIT_BASES   += onBaitBases;
            this.metrics.NEAR_BAIT_BASES += (mappedBases - onBaitBases);            
        }
        else {
            this.metrics.OFF_BAIT_BASES += mappedBases;
        }

    }

    /** Calculates a few last summary metrics and then returns the metrics calculated. */
    public HsMetrics getMetrics() {
        this.metrics.calculateDerivedMetrics();
        calculateTargetCoverageMetrics();
        return this.metrics;
    }

    /** Calculates how much additional sequencing is needed to raise 80% of bases to the mean for the lane. */
    private void calculateTargetCoverageMetrics() {
        short[] depths = new short[(int) this.metrics.TARGET_TERRITORY];  // may not use entire array
        int zeroCoverageTargets = 0;
        int depthIndex = 0;
        double totalCoverage = 0;
        int basesConsidered = 0;

        for (Coverage c : this.coverageByTarget.values()) {
            if (!c.hasCoverage()) {
                ++zeroCoverageTargets;
                continue;
            }

            final short[] targetDepths = c.getDepths();
            basesConsidered += targetDepths.length;

            for (short depth : targetDepths) {
                depths[depthIndex++] = depth;
                totalCoverage += depth;
            }
        }

        this.metrics.MEAN_TARGET_COVERAGE = totalCoverage / basesConsidered;

        // Sort the array (ASCENDING) and then find the base the coverage value that lies at the 80%
        // line, which is actually at 20% into the array now
        Arrays.sort(depths);
        int indexOf80thPercentile = (depths.length - basesConsidered) + (int) (basesConsidered * 0.2);
        int coverageAt80thPercentile = depths[indexOf80thPercentile];
        this.metrics.FOLD_80_BASE_PENALTY = this.metrics.MEAN_TARGET_COVERAGE / coverageAt80thPercentile;
        this.metrics.ZERO_CVG_TARGETS_PCT = zeroCoverageTargets / (double) this.targets.getIntervals().size();
    }
}
