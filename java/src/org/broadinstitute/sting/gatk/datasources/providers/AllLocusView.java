package org.broadinstitute.sting.gatk.datasources.providers;

import java.util.NoSuchElementException;
import java.util.ArrayList;

import org.broadinstitute.sting.gatk.iterators.GenomeLocusIterator;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import net.sf.samtools.SAMRecord;
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
 * A LocusView over which the user can iterate.  
 */

public class AllLocusView extends LocusView {
    private GenomeLocusIterator locusIterator;

    /**
     * Gets the current position in the view.
     */
    private GenomeLoc nextPosition = null;

    /**
     * What's the context for the last locus accessed?
     * @param provider
     */
    private LocusContext nextLocusContext = null;

    /**
     * Create a new queue of locus contexts.
     * @param provider
     */
    public AllLocusView(ShardDataProvider provider) {                
        super( provider );
        // Seed the state tracking members with the first possible seek position and the first possible locus context.
        locusIterator = new GenomeLocusIterator( provider.getShard().getGenomeLoc() );
        if( locusIterator.hasNext() ) {
            nextPosition = locusIterator.next();
            nextLocusContext = hasNextLocusContext() ? nextLocusContext() : createEmptyLocusContext(nextPosition);
        }
    }

    public boolean hasNext() {
        return nextPosition != null;
    }

    public LocusContext next() {
        GenomeLoc currentPosition = nextPosition;
        if( currentPosition == null )
            throw new NoSuchElementException("No next is available in the all locus view");

        // Determine the next locus.
        nextPosition = locusIterator.hasNext() ? locusIterator.next() : null;

        // Crank the iterator to (if possible) or past the next context.
        while( nextLocusContext != null && nextLocusContext.getLocation().isBefore(currentPosition) && hasNextLocusContext() )
            nextLocusContext = nextLocusContext();

        // If actual data is present, return it.  Otherwise, return empty data.
        if( nextLocusContext != null && nextLocusContext.getLocation().equals(currentPosition) )
            return nextLocusContext;
        else
            return createEmptyLocusContext( currentPosition );
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
