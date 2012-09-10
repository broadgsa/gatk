package org.broadinstitute.sting.utils.nanoScheduler;

/**
 * A function that combines a value of MapType with an existing ReduceValue into a new ResultType
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 9:49 AM
 */
public interface NSReduceFunction<MapType, ReduceType> {
    /**
     * Combine one with sum into a new ReduceType
     * @param one the result of a map call on an input element
     * @param sum the cumulative reduce result over all previous map calls
     * @return
     */
    public ReduceType apply(MapType one, ReduceType sum);
}
