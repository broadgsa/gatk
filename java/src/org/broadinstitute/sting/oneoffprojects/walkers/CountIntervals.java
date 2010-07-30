package org.broadinstitute.sting.oneoffprojects.walkers;

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

/**
 * Counts the number of contiguous regions the walker traverses over. Slower than it needs to be, but
 * very useful since overlapping intervals get merged, so you can count the number of intervals the GATK merges down to.
 * This was its very first use.
 */
public class CountIntervals extends RefWalker<GenomeLoc, Pair<GenomeLoc,Long>> {

    public Pair<GenomeLoc,Long> reduceInit() {
        return new Pair<GenomeLoc,Long>(null,0l);
    }

    public GenomeLoc map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) {
            return null;
        }

        return ref.getLocus();
    }

    public Pair<GenomeLoc,Long> reduce(GenomeLoc loc, Pair<GenomeLoc,Long> prev) {
        if ( prev.first == null || prev.first.distance(loc) > 1 ) {
            prev.second ++;
        }

        prev.first = loc;

        return prev;
    }

    public void onTraversalDone(Pair<GenomeLoc,Long> finalReduce ) {
        out.printf("Number of contiguous intervals: %d",finalReduce.second);
    }
}
