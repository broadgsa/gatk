package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.ArrayList;
import java.util.List;

/**
 * Same as CountRod but this one is a reference walker
 */
public class CountRodByRefWalker extends RefWalker<CountRodWalker.Datum, Pair<ExpandingArrayList<Long>, Long>> {
    @Argument(fullName = "verbose", shortName = "v", doc="If true, Countrod will print out detailed information about the rods it finds and locations", required = false)
    public boolean verbose = false;

    @Argument(fullName = "showSkipped", shortName = "s", doc="If true, CountRod will print out the skippped locations", required = false)
    public boolean showSkipped = false;

    CountRodWalker crw = new CountRodWalker();

    public void initialize() {
        crw.verbose = verbose;
        crw.showSkipped = showSkipped;
        crw.out = out;
        crw.err = err;
    }

    public CountRodWalker.Datum map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return crw.map(tracker, ref, context);
    }

    public Pair<ExpandingArrayList<Long>, Long> reduceInit() {
        return crw.reduceInit();
    }

    public Pair<ExpandingArrayList<Long>, Long> reduce(CountRodWalker.Datum point, Pair<ExpandingArrayList<Long>, Long> sum) {
        return crw.reduce(point, sum);
    }
}