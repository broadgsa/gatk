/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam.util;

import java.util.Iterator;

/**
 * PeekIterator is a better class to use than this.
 * @param <T>
 * @param <ITERATOR>
 */
public class NonDestructiveIterator<T, ITERATOR extends Iterator<T>> {
    private T current = null;
    private final ITERATOR underlyingIterator;

    public NonDestructiveIterator(final ITERATOR underlyingIterator) {
        this.underlyingIterator = underlyingIterator;
        advance();
    }

    public T getCurrent() {
        return current;
    }

    public ITERATOR getUnderlyingIterator() {
        return underlyingIterator;
    }

    public boolean advance() {
        if (this.underlyingIterator.hasNext()) {
            current = this.underlyingIterator.next();
        } else {
            current = null;
        }
        return hasCurrent();
    }

    public boolean hasCurrent() {
        return getCurrent() != null;
    }
}
