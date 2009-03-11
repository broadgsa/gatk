package org.broadinstitute.sting.utils;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Feb 24, 2009
 * Time: 10:15:19 AM
 * To change this template use File | Settings | File Templates.
 */
public interface Predicate<T> {
    public boolean apply(T arg);
}

