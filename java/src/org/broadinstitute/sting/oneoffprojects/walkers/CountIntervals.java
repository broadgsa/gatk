package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;

import java.util.List;
import java.io.PrintStream;

/**
 * Counts the number of contiguous regions the walker traverses over. Slower than it needs to be, but
 * very useful since overlapping intervals get merged, so you can count the number of intervals the GATK merges down to.
 * This was its very first use.
 */
public class CountIntervals extends RefWalker<Long, Long> {
    @Output
    PrintStream out;

    @Argument(fullName="numOverlaps",shortName="no",doc="Count all occurrences of X or more overlapping intervals; defaults to 2", required=false)
    int numOverlaps = 2;

    public Long reduceInit() {
        return 0l;
    }

    public boolean isReduceByInterval() { return true; }

    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) {
            return null;
        }

        List<GATKFeature> checkIntervals = tracker.getGATKFeatureMetaData("check",false);
        return (long) checkIntervals.size();
    }

    public Long reduce(Long loc, Long prev) {
        if ( loc == null ) {
            return 0l;
        } else {
            return Math.max(prev,loc);
        }
    }

    public void onTraversalDone(List<Pair<GenomeLoc,Long>> finalReduce) {
        long count = 0;
        for ( Pair<GenomeLoc,Long> g : finalReduce ) {
            if ( g.second >= numOverlaps) {
                count ++;
            }
        }
        out.printf("Number of contiguous intervals: %d",count);
    }
}
