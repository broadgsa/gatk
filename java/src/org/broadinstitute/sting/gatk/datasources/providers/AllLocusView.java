package org.broadinstitute.sting.gatk.datasources.providers;

import java.util.NoSuchElementException;
import java.util.ArrayList;
import java.util.Collections;

import org.broadinstitute.sting.gatk.iterators.GenomeLocusIterator;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
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
     * Gets the next position in the view: next call to next() will jump there.
     * Note that both nextPosition and nextLocus are PRE-read and cached.
     */
    private GenomeLoc nextPosition = null;

    /**
     * What's the next available context?
     */
    private AlignmentContext nextLocus = null;

    /**
     * Create a new queue of locus contexts.
     * @param provider
     */
    public AllLocusView(LocusShardDataProvider provider) {                
        super( provider );
        // Seed the state tracking members with the first possible seek position and the first possible locus context.
        locusIterator = new GenomeLocusIterator(genomeLocParser,provider.getLocus());
        if( locusIterator.hasNext() ) {
            // cache next position and next alignment context
            nextPosition = locusIterator.next();
            nextLocus = hasNextLocus() ? nextLocus() : createEmptyLocus(nextPosition);
        }
    }

    public boolean hasNext() {
        return nextPosition != null;
    }

    public AlignmentContext next() {

        GenomeLoc currentPosition = nextPosition;
        if( currentPosition == null )
            throw new NoSuchElementException("No next is available in the all locus view");

        // Crank the iterator to (if possible) or past the next context.
        while( nextLocus != null && nextLocus.getLocation().isBefore(currentPosition) && hasNextLocus() )
            nextLocus = nextLocus();

        AlignmentContext currentLocus = null; // context we are gonna return

        // If actual data is present, return it.  Otherwise, return empty data.
        if( nextLocus != null && nextLocus.getLocation().equals(currentPosition) ) {
            currentLocus = nextLocus; // found alignment context at the current position
            nextLocus = hasNextLocus() ? nextLocus() : null; 
        }
        else
            currentLocus = createEmptyLocus( currentPosition );

        // Determine the next locus. The trick is that we may have more than one alignment context at the same
        // reference position (regular base pileup, then extended pileup). If next alignment context (that we just pre-read)
        // is still at the current position, we do not increment current position and wait for next call to next() to return
        // that context. If we know that next context is past the current position, we are done with current
        // position
        if ( nextLocus == null || ! nextLocus.getLocation().equals(currentPosition) )
            nextPosition = locusIterator.hasNext() ? locusIterator.next() : null;


        return currentLocus;
    }

    /**
     * Creates a blank locus context at the specified location.
     * @param site Site at which to create the blank locus context.
     * @return empty context.
     */
    private AlignmentContext createEmptyLocus( GenomeLoc site ) {        
        return new AlignmentContext(site,new ReadBackedPileupImpl(site,new ArrayList<SAMRecord>(), new ArrayList<Integer>()));
    }
}
