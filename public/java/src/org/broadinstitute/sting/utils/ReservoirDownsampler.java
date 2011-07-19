package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

/**
 * Randomly downsample from a stream of elements.  This algorithm is a direct,
 * naive implementation of reservoir downsampling as described in "Random Downsampling
 * with a Reservoir" (Vitter 1985).  At time of writing, this paper is located here:
 * http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.138.784&rep=rep1&type=pdf

 * @author mhanna
 * @version 0.1
 */
public class ReservoirDownsampler<T> {
    /**
     * The reservoir of elements tracked by this downsampler.
     */
    private final ArrayList<T> reservoir;

    /**
     * What is the maximum number of reads that can be returned in a single batch.
     */
    private final int maxElements;

    /**
     * Create a new downsampler with the given source iterator and given comparator.
     * @param maxElements What is the maximum number of reads that can be returned in any call of this
     */
    public ReservoirDownsampler(final int maxElements) {
        if(maxElements < 0)
            throw new ReviewedStingException("Unable to work with an negative size collection of elements");
        this.reservoir = new ArrayList<T>(maxElements);
        this.maxElements = maxElements;
    }

    /**
     * Returns the eliminated element.
     * @param element Eliminated element; null if no element has been eliminated.
     * @return
     */
    public T add(T element) {
        if(maxElements <= 0)
            return element;
        else if(reservoir.size() < maxElements) {
            reservoir.add(element);
            return null;
        }
        else {
            // Get a uniformly distributed int. If the chosen slot lives within the partition, replace the entry in that slot with the newest entry.
            int slot = GenomeAnalysisEngine.getRandomGenerator().nextInt(maxElements);
            if(slot >= 0 && slot < maxElements) {
                T displaced = reservoir.get(slot);
                reservoir.set(slot,element);
                return displaced;
            }
            else
                return element;
        }
    }

    public boolean addAll(Collection<? extends T> elements) {
        boolean added = false;
        for(T element: elements)
            added |= (add(element) != null);
        return added;
    }

    /**
     * Returns the contents of this reservoir, downsampled to the given value.  Note that the return value
     * @return The downsampled contents of this reservoir.
     */
    public Collection<T> getDownsampledContents() {
        return reservoir;
    }

    public void clear() {
        reservoir.clear();
    }

    public boolean isEmpty() {
        return reservoir.isEmpty();
    }

    public int size() {
        return reservoir.size();
    }

    public Iterator<T> iterator() {
        return reservoir.iterator();
    }

    public boolean contains(Object o) {
        return reservoir.contains(o);
    }

    public boolean containsAll(Collection<?> elements) {
        return reservoir.containsAll(elements);
    }

    public boolean retainAll(Collection<?> elements) {
        return reservoir.retainAll(elements);
    }

    public boolean remove(Object o) {
        return reservoir.remove(o);
    }

    public boolean removeAll(Collection<?> elements) {
        return reservoir.removeAll(elements);
    }

    public Object[] toArray() {
        Object[] contents = new Object[reservoir.size()];
        reservoir.toArray(contents);
        return contents;
    }

    public <T> T[] toArray(T[] array) {
        return reservoir.toArray(array);
    }
}
