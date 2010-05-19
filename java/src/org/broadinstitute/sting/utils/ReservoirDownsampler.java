package org.broadinstitute.sting.utils;

import net.sf.picard.util.PeekableIterator;

import java.util.*;

/**
 * Randomly downsample from a stream of elements.  This algorithm is a direct,
 * naive implementation of reservoir downsampling as described in "Random Downsampling
 * with a Reservoir" (Vitter 1985).  At time of writing, this paper is located here:
 * http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.138.784&rep=rep1&type=pdf

 * @author mhanna
 * @version 0.1
 */
public class ReservoirDownsampler<T> implements Collection<T> {
    /**
     * Create a random number generator with a random, but reproducible, seed.
     */
    private final Random random = new Random(47382911L);

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
            throw new StingException("Unable to work with an negative size collection of elements");
        this.reservoir = new ArrayList<T>(maxElements);
        this.maxElements = maxElements;
    }

    @Override
    public boolean add(T element) {
        if(maxElements <= 0)
            return false;
        else if(reservoir.size() < maxElements) {
            reservoir.add(element);
            return true;
        }
        else {
            // Get a uniformly distributed int. If the chosen slot lives within the partition, replace the entry in that slot with the newest entry.
            int slot = random.nextInt(maxElements);
            if(slot >= 0 && slot < maxElements) {
                reservoir.set(slot,element);
                return true;
            }
            else
                return false;
        }
    }

    @Override
    public boolean addAll(Collection<? extends T> elements) {
        boolean added = false;
        for(T element: elements)
            added |= add(element);
        return added;
    }

    /**
     * Returns the contents of this reservoir, downsampled to the given value.  Note that the return value
     * @return The downsampled contents of this reservoir.
     */
    public Collection<T> getDownsampledContents() {
        return (Collection<T>)reservoir.clone();
    }

    @Override
    public void clear() {
        reservoir.clear();
    }

    @Override
    public boolean isEmpty() {
        return reservoir.isEmpty();
    }

    @Override
    public int size() {
        return reservoir.size();
    }

    @Override
    public Iterator<T> iterator() {
        return reservoir.iterator();
    }

    @Override
    public boolean contains(Object o) {
        return reservoir.contains(o);
    }

    @Override
    public boolean containsAll(Collection<?> elements) {
        return reservoir.containsAll(elements);
    }

    @Override
    public boolean retainAll(Collection<?> elements) {
        return reservoir.retainAll(elements);
    }

    @Override
    public boolean remove(Object o) {
        return reservoir.remove(o);
    }

    @Override
    public boolean removeAll(Collection<?> elements) {
        return reservoir.removeAll(elements);
    }

    @Override
    public Object[] toArray() {
        Object[] contents = new Object[reservoir.size()];
        reservoir.toArray(contents);
        return contents;
    }

    @Override
    public <T> T[] toArray(T[] array) {
        return reservoir.toArray(array);
    }
}
