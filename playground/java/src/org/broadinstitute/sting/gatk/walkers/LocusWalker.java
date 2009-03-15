package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.LocusContext;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
public interface LocusWalker<MapType, ReduceType> {
    void initialize();
    public String walkerType();

    // Do we actually want to operate on the context?
    boolean filter(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context);

    // Map over the org.broadinstitute.sting.gatk.LocusContext
    MapType map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context);

    // Given result of map function
    ReduceType reduceInit();
    ReduceType reduce(MapType value, ReduceType sum);

    void onTraversalDone();
}
