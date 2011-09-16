package org.broadinstitute.sting.gatk.walkers.qc;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;

import java.io.PrintStream;
import java.util.Collections;
import java.util.List;

/**
 * Counts the number of contiguous regions the walker traverses over. Slower than it needs to be, but
 * very useful since overlapping intervals get merged, so you can count the number of intervals the GATK merges down to.
 * This was its very first use.
 */
public class CountIntervals extends RefWalker<Long, Long> {
    @Output
    PrintStream out;

    @Input(fullName="check", shortName = "check", doc="Any number of RODs", required=false)
    public List<RodBinding<Feature>> features = Collections.emptyList();

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

        List<Feature> checkIntervals = tracker.getValues(features);
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
