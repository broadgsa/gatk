package org.broadinstitute.sting.atk.modules;

import org.broadinstitute.sting.atk.LocusWalker;
import org.broadinstitute.sting.atk.LocusIterator;
import org.broadinstitute.sting.atk.LocusContext;
import org.broadinstitute.sting.utils.ReferenceOrderedDatum;
import net.sf.samtools.SAMRecord;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class BasicLociWalker<MapType, ReduceType> implements LocusWalker<MapType, ReduceType> {
    public void initialize() {
        ;
    }

    public String walkerType() { return "ByLocus"; }

    // Do we actually want to operate on the context?
    public boolean filter(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
        return true;    // We are keeping all the reads
    }

    public void onTraveralDone() {
    }

    // These three capabilities must be overidden
    public abstract MapType map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context);
    public abstract ReduceType reduceInit();
    public abstract ReduceType reduce(MapType value, ReduceType sum);

}
