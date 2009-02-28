/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam.util;

import java.util.Iterator;

public class PeekIterator<T> implements Iterator<T> {
    Iterator<T> underlyingIterator;
    T peekedElement = null;

    public PeekIterator(final Iterator<T> underlyingIterator) {
        this.underlyingIterator = underlyingIterator;
    }

    public boolean hasNext() {
        return peekedElement != null || underlyingIterator.hasNext();  
    }

    public T next() {
        if (peekedElement != null) {
            final T ret = peekedElement;
            peekedElement = null;
            return ret;
        }
        return underlyingIterator.next();
    }

    public T peek() {
        if (peekedElement == null) {
            peekedElement = underlyingIterator.next();
        }
        return peekedElement;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    public Iterator<T> getUnderlyingIterator() {
        return underlyingIterator;
    }
}
