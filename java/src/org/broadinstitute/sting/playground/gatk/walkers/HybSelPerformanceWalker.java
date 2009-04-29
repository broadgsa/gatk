package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;

import net.sf.samtools.SAMRecord;

public class HybSelPerformanceWalker extends LocusWalker<Integer, HybSelPerformanceWalker.TargetInfo> {
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

            // TODO: is there a better way to do this?
            if (read.getNotPrimaryAlignmentFlag() ||
                read.getDuplicateReadFlag() ||
                read.getReadUnmappedFlag() ||
                read.getMappingQuality() <= -1
                    ) {
                continue;
            }
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
        out.println("location\tlength\tavg_coverage\tnormalized_coverage\thit_twice");

        // first zip through and calculate the total average coverage
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

            out.printf("%s:%d-%d\t%d\t%6.4f\t%6.4f\t%d\n",
                       target.getContig(), target.getStart(), target.getStop(), length,
                        avgCoverage, normCoverage, ((ti.hitTwice)?1:0)
                    );


        }
    }
}