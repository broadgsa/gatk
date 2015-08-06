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

package org.broadinstitute.gatk.utils.refdata.utils;

import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Adapter to allow checking if the wrapped iterator was closed.
 * Creating an CCTI also adds it to the list returned from getThreadIterators().
 * @param <T> feature
 */
public class CheckableCloseableTribbleIterator<T extends Feature> implements CloseableTribbleIterator<T> {
    private final CloseableTribbleIterator<T> iterator;
    private boolean closed = false;

    private static ThreadLocal<List<CheckableCloseableTribbleIterator<? extends Feature>>> threadIterators =
            new ThreadLocal<List<CheckableCloseableTribbleIterator<? extends Feature>>>() {
                @Override
                protected List<CheckableCloseableTribbleIterator<? extends Feature>> initialValue() {
                    return new ArrayList<CheckableCloseableTribbleIterator<? extends Feature>>();
                }
            };

    public CheckableCloseableTribbleIterator(CloseableTribbleIterator<T> iterator) {
        this.iterator = iterator;
        threadIterators.get().add(this);
    }

    /**
     * Returns the list of iterators created on this thread since the last time clearCreatedIterators() was called.
     * @return the list of iterators created on this thread since the last time clearCreatedIterators() was called.
     */
    public static List<CheckableCloseableTribbleIterator<? extends Feature>> getThreadIterators() {
        return threadIterators.get();
    }

    /**
     * Clears the tracked list of iterators created on this thread.
     */
    public static void clearThreadIterators() {
        threadIterators.get().clear();
    }

    @Override
    public void close() {
        iterator.close();
        this.closed = true;
    }

    /**
     * Returns true if this iterator was properly closed.
     * @return true if this iterator was properly closed.
     */
    public boolean isClosed() {
        return closed;
    }

    @Override public Iterator<T> iterator() { return this; }
    @Override public boolean hasNext() { return iterator.hasNext(); }
    @Override public T next() { return iterator.next(); }
    @Override public void remove() { iterator.remove(); }
}
