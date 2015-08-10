/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.iterators;

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