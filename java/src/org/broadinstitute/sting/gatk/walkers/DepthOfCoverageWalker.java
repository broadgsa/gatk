package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Pair;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class DepthOfCoverageWalker extends LocusWalker<Integer, Pair<Long, Long>> {
    @Argument(fullName="printall",required=false,defaultValue="true")
    public String printAllLoci;  // booleans don't work

    public Integer map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
        if (printAllLoci.equals("true"))
            out.printf("%s: %d%n", context.getLocation(), context.getReads().size() );
        return context.getReads().size();
    }

    public Pair<Long, Long> reduceInit() { return new Pair<Long, Long>((long)0,(long)0); }

    public Pair<Long, Long> reduce(Integer value, Pair<Long, Long> sum) {
        long left = new Long(value.longValue() + sum.getFirst().longValue());
        long right = new Long(sum.getSecond().longValue() + 1);
        return new Pair<Long, Long>(left, right);
    }

    public void onTraversalDone(Pair<Long, Long> result) {
        out.println("Average depth of coverage is: " + ((double)result.getFirst() / (double)result.getSecond()));
    }
}