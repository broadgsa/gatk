package org.broadinstitute.sting.gatk.dataSources.providers;

import net.sf.samtools.SAMRecord;

import java.util.Iterator;
import java.util.NoSuchElementException;

import edu.mit.broad.picard.filter.FilteringIterator;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.gatk.iterators.LocusContextIteratorByHanger;
import org.broadinstitute.sting.gatk.iterators.LocusContextIterator;
import org.broadinstitute.sting.gatk.iterators.LocusIterator;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
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

public class IterableLocusContextQueue implements LocusContextQueue, LocusIterator {
    private Shard shard;
    private LocusContextIterator loci;

    /**
     * What's the context for the last locus accessed?
     * @param provider
     */
    private LocusContext nextLocusContext = null;    

    /**
     * Create a new queue of locus contexts.
     * @param provider
     */
    public IterableLocusContextQueue(ShardDataProvider provider) {
        Iterator<SAMRecord> reads = new FilteringIterator(provider.getReadIterator(), new TraversalEngine.locusStreamFilterFunc());
        this.loci = new LocusContextIteratorByHanger(reads);
        this.shard = provider.getShard();
    }

    /**
     * Is there another locus present in this iterator.
     * @return True if another locus present in this iterator.  Otherwise, false.
     */
    public boolean hasNext() {
        return loci.hasNext();
    }

    /**
     * Retrieves the next element in the queue.
     * @return Next element in the queue.
     */
    public GenomeLoc next() {
        do {
            nextLocusContext = loci.next();
        }
        while( nextLocusContext.getLocation().isBefore(shard.getGenomeLoc()) );
        
        return nextLocusContext.getLocation();
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
        if( nextLocusContext == null )
            throw new NoSuchElementException("No more elements remaining in queue");
        return nextLocusContext;
    }

    /**
     * Seek to the specified position in the contig.
     * @param seekPoint
     */
    public LocusContextQueue seek( GenomeLoc seekPoint ) {
        if( nextLocusContext == null || !seekPoint.equals(nextLocusContext.getLocation()) ) {
            nextLocusContext = null;
            throw new IllegalArgumentException("IterableLocusContextQueue doesn't support seeking and iterator is in the wrong position.");
        }
        return this;
    }

}
