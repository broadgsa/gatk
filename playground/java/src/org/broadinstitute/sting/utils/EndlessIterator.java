package org.broadinstitute.sting.utils;

import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Feb 24, 2009
 * Time: 10:24:38 AM
 * To change this template use File | Settings | File Templates.
 */
public class EndlessIterator<T> implements Iterator<T> {
    private T value;
    
    public EndlessIterator(T value) {
        this.value = value;
    }

    public boolean hasNext() {
        return true;
    }

    public T next() {
        return this.value;
    }

    public void remove () {
        throw new UnsupportedOperationException();
    }
}