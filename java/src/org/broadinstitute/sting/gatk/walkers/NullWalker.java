package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

// Null traversal. For ATK performance measuring.
// j.maguire 3-7-2009

public class NullWalker extends LocusWalker<Integer, Integer> {
    public void initialize() {
    }

    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return true;    // We are keeping all the reads
    }

    // Map over the org.broadinstitute.sting.gatk.contexts.AlignmentContext
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context)
    {
        return 1;
    }

    // Given result of map function
    public Integer reduceInit() 
    {
        return 0; 
    }
    public Integer reduce(Integer value, Integer sum) 
    {
        return 0;
    }

    public void onTraversalDone() {
    }
}
