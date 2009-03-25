package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class DepthOfCoverageWalker extends BasicLociWalker<Integer, Integer> {
    public Integer map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
        System.out.printf("%s: %d%n", context.getLocation(), context.getReads().size() );
        return 1;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}