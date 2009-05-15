package org.broadinstitute.sting.gatk.dataSources.providers;

import org.broadinstitute.sting.gatk.iterators.LocusContextIterator;
import org.broadinstitute.sting.gatk.iterators.LocusContextIteratorByHanger;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.apache.log4j.Logger;
import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Iterator;

import edu.mit.broad.picard.filter.FilteringIterator;
/**
 * User: hanna
 * Date: May 12, 2009
 * Time: 11:24:42 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A queue of locus contexts.  Provides unidirectional seek.  Stripped down
 * implementation of java.util.Queue interface.
 */

public class SeekableLocusContextQueue extends LocusContextQueue {
    /**
     * Gets the position to which the last seek was requested.
     */
    private GenomeLoc seekPoint;

    /**
     * What's the context for the last locus accessed?
     * @param provider
     */
    private LocusContext nextLocusContext = null;

    private static Logger logger = Logger.getLogger(SeekableLocusContextQueue.class);

    /**
     * Create a new queue of locus contexts.
     * @param provider
     */
    public SeekableLocusContextQueue(ShardDataProvider provider) {
        super(provider);

        // Seed the state tracking members with the first possible seek position and the first possible locus context.
        seekPoint = new GenomeLoc(shard.getGenomeLoc().getContigIndex(),shard.getGenomeLoc().getStart());

        if( hasNextLocusContext() )
            nextLocusContext = getNextLocusContext();
        else
            nextLocusContext = this.createEmptyLocusContext(seekPoint);                
    }

    /**
     * Get the locus context at the given position.
     * @return Locus context, or null if no locus context exists at this position.
     */
    public LocusContext peek() {
        // Haven't reached the next locus context in the list yet.  Return null.
        if( seekPoint.isBefore(nextLocusContext.getLocation()) )
            return createEmptyLocusContext(seekPoint);

        return nextLocusContext;
    }

    /**
     * Seek to the given point the queue of locus contexts.
     * @param target Target base pair to which to seek.  Must be a single base pair.
     * @return an instance of itself for parameter chaining.
     */
    public LocusContextQueue seek(GenomeLoc target) {
        if( !target.isSingleBP() )
            throw new IllegalArgumentException("Seek point must be a single base pair.");

        // If outside the range of the target, throw an illegal argument exception.
        if( target.isBefore(shard.getGenomeLoc()) || target.isPast(shard.getGenomeLoc()))
            throw new IllegalArgumentException(String.format("Target is out of range; target = %s, valid range = %s",target,shard.getGenomeLoc()));

        seekPoint = (GenomeLoc)target.clone();

        // Search for the next locus context following the target positions.
        while (nextLocusContext.getLocation().isBefore(target) && hasNextLocusContext() ) {
            logger.debug(String.format("  current locus is %s vs %s => %d", nextLocusContext.getLocation(),
                                                                            target,
                                                                            nextLocusContext.getLocation().compareTo(target)));
            nextLocusContext = getNextLocusContext();
        }

        // Couldn't find a next?  Force the nextLocusContext to null.
        if( nextLocusContext.getLocation().isBefore(target) && !hasNextLocusContext() )
            nextLocusContext = createEmptyLocusContext( seekPoint );

        return this;
    }

    /**
     * Gets the point to which the queue has currently seeked.
     * @return Single bp position where the queue has been positioned.  A locus context may or may not
     *         exist at this point.
     */
    public GenomeLoc getSeekPoint() {
        return seekPoint;
    }

    /**
     * Creates a blank locus context at the specified location.
     * @param site Site at which to create the blank locus context.
     * @return empty context.
     */
    private LocusContext createEmptyLocusContext( GenomeLoc site ) {
        return new LocusContext(site, new ArrayList<SAMRecord>(), new ArrayList<Integer>());
    }
}
