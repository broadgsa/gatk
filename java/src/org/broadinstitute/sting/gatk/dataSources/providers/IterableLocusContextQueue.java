package org.broadinstitute.sting.gatk.dataSources.providers;

import java.util.NoSuchElementException;

import org.broadinstitute.sting.gatk.iterators.LocusIterator;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.GenomeLoc;
/**
 * User: hanna
 * Date: May 13, 2009
 * Time: 3:32:30 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A LocusContextQueue over which the user can iterate.  
 */

public class IterableLocusContextQueue extends LocusContextQueue implements LocusIterator {
    /**
     * What's the context for the last locus accessed?
     * @param provider
     */
    private LocusContext prefetched = null;

    /**
     * Has this prefetch been consumed?  If this flag is set,
     * the prefetch will skip to the next argument in the system.
     */
    private boolean prefetchConsumed = true;

    /**
     * Create a new queue of locus contexts.
     * @param provider
     */
    public IterableLocusContextQueue(ShardDataProvider provider) {
        super( provider );
    }

    /**
     * Is there another locus present in this iterator.
     * @return True if another locus present in this iterator.  Otherwise, false.
     */
    public boolean hasNext() {
        prefetchLocusContext();
        return prefetched != null;
    }

    /**
     * Retrieves the next element in the queue.
     * @return Next element in the queue.
     */
    public GenomeLoc next() {
        prefetchLocusContext();
        prefetchConsumed = true;
        // Signal that the prefetcher needs to grab another entry off the queue.        
        return prefetched.getLocation();
    }

    /**
     * Find the next locus context within the bounds of a member variable and store
     * it in the prefetched member variable.  When the prefetch is consumed, the 'consumer'
     * should signal it as such by marking prefetchConsumed = true.
     */
    private void prefetchLocusContext() {
        if( !prefetchConsumed )
            return;

        prefetched = null;
        prefetchConsumed = false;        

        // If another locus context bounded by this shard exists, find it.
        boolean prefetchOutOfBounds = true;
        while( hasNextLocusContext() && prefetchOutOfBounds ) {
            prefetched = getNextLocusContext();
            prefetchOutOfBounds = (prefetched.getLocation().isBefore(shard.getGenomeLoc()) ||
                                   prefetched.getLocation().isPast(shard.getGenomeLoc()));
        }

        // Can't find a valid prefetch?  Set prefetch to null.  If prefetched == null and
        // prefetchConsumed == false, the queue is out of entries.
        if( prefetchOutOfBounds )
            prefetched = null;
    }

    /**
     * Unsupported.
     */
    public void remove() {
        throw new UnsupportedOperationException("Unable to remove elements from this queue.");
    }

    /**
     * Peek at the next locus context in the chain.
     * @return
     */
    public LocusContext peek() {
        if( prefetched == null )
            throw new NoSuchElementException("No more elements remaining in queue");
        return prefetched;
    }

    /**
     * Seek to the specified position in the contig.
     * @param seekPoint
     */
    public LocusContextQueue seek( GenomeLoc seekPoint ) {
        if( prefetched == null || !seekPoint.equals(prefetched.getLocation()) )
            throw new IllegalArgumentException("IterableLocusContextQueue doesn't support seeking and iterator is in the wrong position.");
        return this;
    }

}
