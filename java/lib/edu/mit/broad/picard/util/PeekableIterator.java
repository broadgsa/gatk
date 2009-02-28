/*
  * The Broad Institute
  * SOFTWARE COPYRIGHT NOTICE AGREEMENT
  * This software and its documentation are copyright Jan 22, 2009 by the
  * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
  *
  * This software is supplied without any warranty or guaranteed support whatsoever. Neither
  * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
  */
package edu.mit.broad.picard.util;

import edu.mit.broad.sam.util.CloseableIterator;

/**
 * Generic Closable Iterator that allows you to peek at the next value before calling next
 */
public class PeekableIterator<Object> implements CloseableIterator<Object> {
    private CloseableIterator<Object> iterator;
    private Object nextObject;

    /** Constructs a new iterator that wraps the supplied iterator. */
    public PeekableIterator(CloseableIterator<Object> iterator) {
        this.iterator = iterator;
        advance();
    }

    /** Closes the underlying iterator. */
    public void close() {
        this.iterator.close();
    }

    /** True if there are more items, in which case both next() and peek() will return a value. */
    public boolean hasNext() {
        return this.nextObject != null;
    }

    /** Returns the next object and advances the iterator. */
    public Object next() {
        Object retval = this.nextObject;
        advance();
        return retval;
    }

    /**
     * Returns the next object but does not advance the iterator. Subsequent calls to peek()
     * and next() will return the same object.
     */
    public Object peek(){
        return this.nextObject;
    }

    private void advance(){
        if (this.iterator.hasNext()) {
            this.nextObject = iterator.next();
        }
        else {
            this.nextObject = null;
        }
    }

    /** Unsupported Operation. */
    public void remove() {
        throw new UnsupportedOperationException("Not supported: remove");
    }
}
