package org.broadinstitute.sting.utils.nanoScheduler;

/**
 * A function that maps from InputType -> ResultType
 *
 * For use with the NanoScheduler
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 9:49 AM
 */
public interface NSMapFunction<InputType, ResultType> {
    /**
     * Return function on input, returning a value of ResultType
     * @param input
     * @return
     */
    public ResultType apply(final InputType input);
}
