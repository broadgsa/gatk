/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sting.utils;

import java.util.Iterator;

public class PushbackIterator<T> implements Iterator<T> {
    Iterator<T> underlyingIterator;
    T pushedElement = null;

    public PushbackIterator(final Iterator<T> underlyingIterator) {
        this.underlyingIterator = underlyingIterator;
    }

    public boolean hasNext() {
        return pushedElement != null || underlyingIterator.hasNext();
    }

    public T next() {
        if (pushedElement != null) {
            final T ret = pushedElement;
            pushedElement = null;
            return ret;
        }
        return underlyingIterator.next();
    }

    public void pushback(T elt) {
        pushedElement = elt;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    public Iterator<T> getUnderlyingIterator() {
        return underlyingIterator;
    }
}