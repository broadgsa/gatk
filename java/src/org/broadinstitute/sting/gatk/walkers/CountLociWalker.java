package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class CountLociWalker extends LocusWalker<Integer, Integer> implements TreeReducible<Integer> {
    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        return 1;
    }

    public Integer reduceInit() { return 0; }

    /**
     * Reduces two subtrees together.  In this case, the implementation of the tree reduce
     * is exactly the same as the implementation of the single reduce.
     */
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}
