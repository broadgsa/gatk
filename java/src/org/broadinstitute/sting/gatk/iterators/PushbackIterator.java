/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package org.broadinstitute.sting.gatk.iterators;

import java.util.Iterator;

public class PushbackIterator<T> implements Iterator<T>, Iterable<T> {
    Iterator<T> underlyingIterator;
    T pushedElement = null;

    public PushbackIterator(final Iterator<T> underlyingIterator) {
        this.underlyingIterator = underlyingIterator;
    }

    public boolean hasNext() {
        return pushedElement != null || underlyingIterator.hasNext();
    }

    public Iterator<T> iterator() {
        return this;
    }

    /**
     * Retrieves, but does not remove, the head of this iterator.
     * @return T the next element in the iterator
     */
    public T element() {
        T x = next();
        pushback(x);
        return x;
    }

    /**
     * @return the next element in the iteration.
     */
    public T next() {
        if (pushedElement != null) {
            final T ret = pushedElement;
            pushedElement = null;
            return ret;
        } else {
            return underlyingIterator.next();
        }
    }

    public void pushback(T elt) {
        assert(pushedElement == null);
        
        pushedElement = elt;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    public Iterator<T> getUnderlyingIterator() {
        return underlyingIterator;
    }
}