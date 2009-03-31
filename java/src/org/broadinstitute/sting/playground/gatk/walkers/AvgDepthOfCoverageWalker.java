package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import java.util.List;

@WalkerName("Average_Depth_Of_Coverage")
public class AvgDepthOfCoverageWalker extends LocusWalker<Integer, Pair<Long, Long>> {

    public Integer map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
        return context.getReads().size();
    }

    public Pair<Long, Long> reduceInit() { return new Pair<Long, Long>((long)0,(long)0); }

    public Pair<Long, Long> reduce(Integer value, Pair<Long, Long> sum) {
        long left = new Long(value.longValue() + sum.getFirst().longValue());
        long right = new Long(sum.getSecond().longValue() + 1);
        return new Pair<Long, Long>(left, right);
    }

    public void onTraversalDone(Pair<Long, Long> result) {
        System.out.println("Average depth of coverage is: " + ((double)result.getFirst() / (double)result.getSecond()));
    }
}
