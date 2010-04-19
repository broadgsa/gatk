package org.broadinstitute.sting.utils;

import net.sf.picard.util.PeekableIterator;

import java.util.*;

/**
 * Randomly downsample from a stream of elements.  This algorithm is a direct,
 * naive implementation of reservoir downsampling as described in "Random Downsampling
 * with a Reservoir" (Vitter 1985).  At time of writing, this paper is located here:
 * http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.138.784&rep=rep1&type=pdf
 *
 * Note that using the ReservoirDownsampler will leave the given iterator in an undefined
 * state.  Do not attempt to use the iterator (other than closing it) after the Downsampler
 * completes.
 *
 * @author mhanna
 * @version 0.1
 */
public class ReservoirDownsampler<T> implements Iterator<Collection<T>> {
    /**
     * Create a random number generator with a random, but reproducible, seed.
     */
    private final Random random = new Random(47382911L);

    /**
     * The data source, wrapped in a peekable input stream.
     */
    private final PeekableIterator<T> iterator;

    /**
     * Used to identify whether two elements are 'equal' in the eyes of the downsampler.
     */
    private final Comparator<T> comparator;

    /**
     * What is the maximum number of reads that can be returned in a single batch.
     */
    private final int maxElements;

    /**
     * Create a new downsampler with the given source iterator and given comparator.
     * @param iterator Source of the data stream.
     * @param comparator Used to compare two records to see whether they're 'equal' at this position.
     * @param maxElements What is the maximum number of reads that can be returned in any call of this
     */
    public ReservoirDownsampler(final Iterator<T> iterator, final Comparator<T> comparator, final int maxElements) {
        this.iterator = new PeekableIterator<T>(iterator);
        this.comparator = comparator;
        if(maxElements < 0)
            throw new StingException("Unable to work with an negative size collection of elements");
        this.maxElements = maxElements;
    }

    public boolean hasNext() {
        return iterator.hasNext();
    }

    /**
     * Gets a collection of 'equal' elements, as judged by the comparator.  If the number of equal elements
     * is greater than the maximum, then the elements in the collection should be a truly random sampling.
     * @return Collection of equal elements.
     */
    public Collection<T> next() {
        if(!hasNext())
            throw new NoSuchElementException("No next element is present.");

        List<T> batch = new ArrayList<T>(maxElements);
        int currentElement = 0;

        // Determine our basis of equality.
        T first = iterator.next();
        if(maxElements > 0)
            batch.add(first);
        currentElement++;

        // Fill the reservoir
        while(iterator.hasNext() &&
              currentElement < maxElements &&
              comparator.compare(first,iterator.peek()) == 0) {
            batch.add(iterator.next());
            currentElement++;
        }

        // Trim off remaining elements, randomly selecting them using the process as described by Vitter.
        while(iterator.hasNext() && comparator.compare(first,iterator.peek()) == 0) {
            T candidate = iterator.next();
            final int slot = random.nextInt(currentElement);
            if(slot >= 0 && slot < maxElements)
                batch.set(slot,candidate);
            currentElement++;
        }

        return batch;
    }

    /**
     * Unsupported; throws exception to that effect.
     */
    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a ReservoirDownsampler.");
    }
}
