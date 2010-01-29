package org.broadinstitute.sting.gatk.iterators;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.NoSuchElementException;
import java.util.Iterator;
import java.util.List;
/**
 * User: hanna
 * Date: May 12, 2009
 * Time: 10:52:47 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Iterates through all of the loci provided in the reference.
 */
public class GenomeLocusIterator implements Iterator<GenomeLoc> {
    /**
     * An iterator to the entire data structure over which we're iterating.
     */
    private final Iterator<GenomeLoc> locusIterator;

    /**
     * The multi-base pair long locus referring to the current locus.
     */
    private GenomeLoc currentLocus = null;

    /**
     * The 1 base pair long location.
     */
    private GenomeLoc currentLocation = null;

    /**
     * Creates an iterator that can traverse over the entire
     * reference specified in the given ShardDataProvider.
     * @param loci the list of loci over which to iterate.
     */
    public GenomeLocusIterator( List<GenomeLoc> loci ) {
        this.locusIterator = loci.iterator();
        seedNextLocus();
    }

    /**
     * Is the iterator still within the locus?
     * @return True if the iterator has more elements.  False otherwise. 
     */
    public boolean hasNext() {
        return currentLocation != null;
    }

    /**
     * Get the next single-base locus context bounded by the iterator.
     * @return GenomeLoc representing the next single-base locus context.
     */
    public GenomeLoc next() {
        if( !hasNext() )
            throw new NoSuchElementException("No elements remaining in bounded reference region.");
        GenomeLoc toReturn = currentLocation.clone();
        seedNextLocus();
        return toReturn;
    }

    public void remove() {
        throw new UnsupportedOperationException( "ReferenceLocusIterator is read-only" );
    }

    /**
     * Position currentLocation at the next locus, if possible.
     */
    private void seedNextLocus() {
        if(currentLocus != null && currentLocation != null)
            currentLocation = GenomeLocParser.incPos(currentLocation);

        // If initializing or the location was pushed off the current locus, reinitialize using the next locus.
        if(currentLocus == null || currentLocation == null || currentLocation.isPast(currentLocus)) {
            currentLocus = currentLocation = null;
            if(locusIterator.hasNext()){
                currentLocus = locusIterator.next();
                currentLocation = GenomeLocParser.createGenomeLoc(currentLocus.getContig(),currentLocus.getStart());
            }
        }
    }
}
