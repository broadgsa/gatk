package org.broadinstitute.sting.utils.nanoScheduler;

/**
 * A function that maps from InputType -> ResultType
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 9:49 AM
 */
public interface ReduceFunction<MapType, ReduceType> {
    public ReduceType init();
    public ReduceType apply(MapType one, ReduceType sum);
}
