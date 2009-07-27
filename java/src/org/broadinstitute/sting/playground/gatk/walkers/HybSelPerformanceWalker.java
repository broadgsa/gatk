package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.util.List;
import java.util.Collection;
import java.io.IOException;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import net.sf.picard.reference.ReferenceSequence;
import edu.mit.broad.picard.util.Interval;
import edu.mit.broad.picard.util.OverlapDetector;

@By(DataSource.REFERENCE)
public class HybSelPerformanceWalker extends LocusWalker<Integer, HybSelPerformanceWalker.TargetInfo> {
    @Argument(fullName="min_mapq", shortName="mmq", required=false, doc="Minimum mapping quality of reads to consider")
    public Integer MIN_MAPQ = 1;

    @Argument(fullName="include_duplicates", shortName="idup", required=false, doc="consider duplicate reads")
    public boolean INCLUDE_DUPLICATE_READS = false;

    @Argument(fullName="free_standing_distance", shortName="fsd", required=false, doc="minimum distance to next interval to consider freestanding")
    public Integer FREE_STANDING_DISTANCE = 500;

    public static class TargetInfo {
        public int counts = 0;

        // did at least two reads hit this target
        public boolean hitTwice = false;

        // TODO: track max and min?
        // TODO: median rather than average?
        // TODO: bin into segments? (requires knowing position)
    }

//    @Argument(fullName="suppressLocusPrinting",required=false,defaultValue="false")
//    public boolean suppressPrinting;

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        List<SAMRecord> reads = context.getReads();

        int depth = 0;
        for ( int i = 0; i < reads.size(); i++ )
        {
            SAMRecord read = reads.get(i);

            if (read.getNotPrimaryAlignmentFlag()) { continue; }
            if (read.getReadUnmappedFlag()) { continue; }
            if (!INCLUDE_DUPLICATE_READS && read.getDuplicateReadFlag()) { continue; }
            if (read.getMappingQuality() < MIN_MAPQ) { continue; }
            depth++;
        }

        return depth;
    }

    /**
     * Return true if your walker wants to reduce each interval separately.  Default is false.
     *
     * If you set this flag, several things will happen.
     *
     * The system will invoke reduceInit() once for each interval being processed, starting a fresh reduce
     * Reduce will accumulate normally at each map unit in the interval
     * However, onTraversalDone(reduce) will be called after each interval is processed.
     * The system will call onTraversalDone( GenomeLoc -> reduce ), after all reductions are done,
     *   which is overloaded here to call onTraversalDone(reduce) for each location
     */
    public boolean isReduceByInterval() {
        return true;
    }

    public TargetInfo reduceInit() { return new TargetInfo(); }

    public TargetInfo reduce(Integer value, TargetInfo sum) {
        sum.counts += value;
        if (value >= 2) { sum.hitTwice = true; }
        return sum;
    }

    public void onTraversalDone(TargetInfo result) {
    }

    @Override
    public void onTraversalDone(List<Pair<GenomeLoc, TargetInfo>> results) {
        out.println("location\tlength\tgc\tavg_coverage\tnormalized_coverage\thit_twice\tfreestanding");

        // first zip through and build an overlap detector of all the intervals, so later
        // we can calculate if this interval is free-standing
        OverlapDetector<Interval> od = new OverlapDetector<Interval>(-FREE_STANDING_DISTANCE,0);
        for(Pair<GenomeLoc, TargetInfo> pair : results) {
            GenomeLoc target = pair.getFirst();
            Interval interval = makeInterval(target);
            od.addLhs(interval, interval);
        }


        // now zip through and calculate the total average coverage
        long totalCoverage = 0;
        long basesConsidered = 0;
        for(Pair<GenomeLoc, TargetInfo> pair : results) {
            GenomeLoc target = pair.getFirst();
            TargetInfo ti = pair.getSecond();

            // as long as it was hit twice, count it
            if(ti.hitTwice) {
                long length = target.getStop() - target.getStart() + 1;
                totalCoverage += ti.counts;
                basesConsidered += length;
            }
        }
        double meanTargetCoverage = totalCoverage / basesConsidered;


        for(Pair<GenomeLoc, TargetInfo> pair : results) {
            GenomeLoc target = pair.getFirst();
            TargetInfo ti = pair.getSecond();
            long length = target.getStop() - target.getStart() + 1;

            double avgCoverage = ((double)ti.counts / (double)length);
            double normCoverage = avgCoverage / meanTargetCoverage;

            // calculate gc for the target
            double gc = calculateGC(target);

            // if there is more than one hit on the overlap detector, it's not freestanding
            Collection<Interval> hits = od.getOverlaps(makeInterval(target));
            boolean freestanding = (hits.size() == 1);

            out.printf("%s:%d-%d\t%d\t%6.4f\t%6.4f\t%6.4f\t%d\t%d\n",
                       target.getContig(), target.getStart(), target.getStop(), length, gc,
                        avgCoverage, normCoverage, ((ti.hitTwice)?1:0), ((freestanding)?1:0)
                    );


        }
    }

    private Interval makeInterval(GenomeLoc target) {
        return new Interval(target.getContig(), (int) target.getStart(), (int) target.getStop());
    }

    private double calculateGC(GenomeLoc target) {
        try {
            IndexedFastaSequenceFile seqFile = new IndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
            ReferenceSequence refSeq = seqFile.getSubsequenceAt(target.getContig(),target.getStart(), target.getStop());


            int gcCount = 0;
            for(char base : StringUtil.bytesToString(refSeq.getBases()).toCharArray()) {
                if (base == 'C' || base == 'c' || base == 'G' || base == 'g') { gcCount++; }
            }
            return ( (double) gcCount ) / ((double) refSeq.getBases().length);

        } catch (IOException ioe) {
            throw new RuntimeException(ioe);
        }

    }
}