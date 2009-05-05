package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
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
    @Argument(fullName="suppressLocusPrinting",doc="Suppress printing",required=false,defaultValue="false")
    public boolean suppressPrinting;

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        if ( !suppressPrinting )
            out.printf("%s: %d%n", context.getLocation(), context.getReads().size() );
        return context.getReads().size();
    }

    public boolean isReduceByInterval() {
        return true;
    }

    public Pair<Long, Long> reduceInit() { return new Pair(0l,0l); }

    public Pair<Long, Long> reduce(Integer value, Pair<Long, Long> sum) {
        long left = value.longValue() + sum.getFirst();
        long right = sum.getSecond() + 1l;
        return new Pair(left, right);
    }

    public void onTraversalDone(Pair<Long, Long> result) {
        out.printf("Average depth of coverage is: %.2f in %d total coverage over %d sites\n", 
                ((double)result.getFirst() / (double)result.getSecond()), result.getFirst(), result.getSecond());
    }
}