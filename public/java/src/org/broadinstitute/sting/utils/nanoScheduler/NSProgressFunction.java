package org.broadinstitute.sting.utils.nanoScheduler;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 9/4/12
 * Time: 2:10 PM
 * To change this template use File | Settings | File Templates.
 */
public interface NSProgressFunction<InputType> {
    public void progress(final InputType lastMapInput);
}
