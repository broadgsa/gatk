package org.broadinstitute.sting.gatk.iterators;

import java.util.Iterator;

public class IterableIterator<T> implements Iterable<T> {
    private Iterator<T> iter;

    public IterableIterator(Iterator<T> iter) {
        this.iter = iter;
    }

    public Iterator<T> iterator() {
        return iter;
    }
}
