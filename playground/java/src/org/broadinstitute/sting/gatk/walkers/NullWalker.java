package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;

import java.util.List;

// Null traversal. For ATK performance measuring.
// j.maguire 3-7-2009

public class NullWalker implements LocusWalker<Integer, Integer> {
    public void initialize() {
    }

    public String walkerType() { return "ByLocus"; }

    // Do we actually want to operate on the context?
    public boolean filter(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
        return true;    // We are keeping all the reads
    }

    // Map over the org.broadinstitute.sting.gatk.LocusContext
    public Integer map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) 
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
