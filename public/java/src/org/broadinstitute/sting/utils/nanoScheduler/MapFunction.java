package org.broadinstitute.sting.utils.nanoScheduler;

/**
 * A function that maps from InputType -> ResultType
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 9:49 AM
 */
public interface MapFunction<InputType, ResultType> {
    public ResultType apply(final InputType input);
}
