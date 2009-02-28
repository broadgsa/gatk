package edu.mit.broad.sting.atk;

import edu.mit.broad.sting.atk.LocusIterator;
import edu.mit.broad.sting.utils.ReferenceOrderedDatum;

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
    boolean filter(List<ReferenceOrderedDatum> rodData, char ref, LocusIterator context);

    // Map over the edu.mit.broad.sting.atk.LocusContext
    MapType map(List<ReferenceOrderedDatum> rodData, char ref, LocusIterator context);

    // Given result of map function
    ReduceType reduceInit();
    ReduceType reduce(MapType value, ReduceType sum);

    void onTraveralDone();
}