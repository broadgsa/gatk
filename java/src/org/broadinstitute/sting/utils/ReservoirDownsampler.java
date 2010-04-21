package org.broadinstitute.sting.utils;

import net.sf.picard.util.PeekableIterator;

import java.util.*;

/**
 * Randomly downsample from a stream of elements.  This algorithm is a direct,
 * naive implementation of reservoir downsampling as described in "Random Downsampling
 * with a Reservoir" (Vitter 1985).  At time of writing, this paper is located here:
 * http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.138.784&rep=rep1&type=pdf
 *
 * Contains an enhancement allowing users to partition downsampled data.  If a partitioner
 * is used, each partition will be allowed to contain maxElements elements.
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
     * Partitions the elements into subsets, each having an equal number of maxElements.
     */
    private final Partitioner<T> partitioner;

    /**
     * What is the maximum number of reads that can be returned in a single batch.
     */
    private final int maxElements;

    /**
     * Create a new downsampler with the given source iterator and given comparator.
     * @param iterator Source of the data stream.
     * @param comparator Used to compare two records to see whether they're 'equal' at this position.
     * @param maxElements What is the maximum number of reads that can be returned in any partition of any call of this iterator.
     */
    public ReservoirDownsampler(final Iterator<T> iterator, final Comparator<T> comparator, final int maxElements) {
        this(iterator,comparator,null,maxElements);
    }

    /**
     * Create a new downsampler with the given source iterator and given comparator.
     * @param iterator Source of the data stream.
     * @param comparator Used to compare two records to see whether they're 'equal' at this position.
     * @param partitioner Used to divide the elements into bins.  Each bin can have maxElements elements.
     * @param maxElements What is the maximum number of reads that can be returned in any call of this
     */
    public ReservoirDownsampler(final Iterator<T> iterator, final Comparator<T> comparator, final Partitioner<T> partitioner, final int maxElements) {
        this.iterator = new PeekableIterator<T>(iterator);
        this.comparator = comparator;
        this.partitioner = partitioner;
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

        Map<Object,Partition<T>> partitions = new HashMap<Object,Partition<T>>();

        // Determine our basis of equality.
        T first = iterator.next();

        if(maxElements > 0)
            getPartitionForEntry(partitions,first).add(first);

        while(iterator.hasNext() && comparator.compare(first,iterator.peek()) == 0) {
            T candidate = iterator.next();
            getPartitionForEntry(partitions,candidate).add(candidate);
        }        

        LinkedList<T> batch = new LinkedList<T>();
        for(Partition<T> partition: partitions.values())
            batch.addAll(partition.elements);

        return batch;
    }

    /**
     * Gets the appropriate partition for the given entry from storage.
     * @param partitions List of partitions from which to choose.
     * @param entry Entry for which to compute the partition.
     * @return The partition associated with this entry.  Will be created if not present.
     */
    private Partition<T> getPartitionForEntry(final Map<Object,Partition<T>> partitions, final T entry) {
        Object partition = partitioner!=null ? partitioner.partition(entry) : null;
        if(!partitions.containsKey(partition))
            partitions.put(partition,new Partition<T>(maxElements));
        return partitions.get(partition);
    }

    /**
     * Unsupported; throws exception to that effect.
     */
    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a ReservoirDownsampler.");
    }

    /**
     * A common interface for a functor that can take data of
     * some type and return an object that can be used to partition
     * that data in some way.  Really just a declaration of a
     * specialized map function.
     */
    public interface Partitioner<T> {
        public Object partition(T input);
    }

    /**
     * Models a partition of a given set of elements.  Knows how to select
     * random elements with replacement.
     * @param <T> Data type for the elements of the partition.
     */
    private class Partition<T> {
        /**
         * How large can this partition grow?
         */
        private final int partitionSize;

        /**
         * The elements of the partition.
         */
        private List<T> elements = new ArrayList<T>();

        /**
         * The total number of elements seen.
         */
        private long elementsSeen = 0;

        public Partition(final int partitionSize) {
            this.partitionSize = partitionSize;
        }

        /**
         * Add a new element to this collection, downsampling as necessary so that the partition
         * stays under partitionSize elements.
         * @param element Element to conditionally add.
         */
        public void add(T element) {
            if(elements.size() < partitionSize)
                elements.add(element);
            else {
                // Get a uniformly distributed long > 0 and remap it to the range from [0,elementsSeen).
                long slot = random.nextLong();
                while(slot == Long.MIN_VALUE)
                    slot = random.nextLong();
                slot = (long)(((float)Math.abs(slot))/Long.MAX_VALUE * (elementsSeen-1));

                // If the chosen slot lives within the partition, replace the entry in that slot with the newest entry.
                if(slot >= 0 && slot < partitionSize)
                    elements.set((int)slot,element);
            }
            elementsSeen++;
        }
    }
}
