package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
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
public abstract class BasicLociWalker<MapType, ReduceType> implements LocusWalker<MapType, ReduceType> {
    public String getName() {
        // Return name of class, trimming 'Walker' from the end if present.
        // TODO: Duplicate of BasicReadWalker.getName().  Eliminate duplication.
        String className = getClass().getSimpleName();
        if(className.endsWith(Walker.class.getSimpleName()))
            return className.substring(0,className.lastIndexOf(Walker.class.getSimpleName()));
        else
            return className;
    }

    public void initialize() {
        ;
    }

    // Do we actually want to operate on the context?
    public boolean filter(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
        return true;    // We are keeping all the reads
    }

    public void onTraversalDone() {
        ;
    }

    // These three capabilities must be overidden
    public abstract MapType map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context);
    public abstract ReduceType reduceInit();
    public abstract ReduceType reduce(MapType value, ReduceType sum);

}
