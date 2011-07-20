package org.broadinstitute.sting.gatk.iterators;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.Iterator;
import java.util.NoSuchElementException;
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
     * Builds individual loci.
     */
    private GenomeLocParser parser;

    /**
     * The entire region over which we're iterating.
     */
    private GenomeLoc completeLocus;

    /**
     * The current position in the traversal.
     */
    private GenomeLoc currentLocus;

    /**
     * Creates an iterator that can traverse over the entire
     * reference specified in the given ShardDataProvider.
     * @param completeLocus Data provider to use as a backing source.
     *                 Provider must have a reference (hasReference() == true).
     */
    public GenomeLocusIterator( GenomeLocParser parser, GenomeLoc completeLocus ) {
        this.parser = parser;
        this.completeLocus = completeLocus;
        this.currentLocus = parser.createGenomeLoc(completeLocus.getContig(),completeLocus.getStart());
    }

    /**
     * Is the iterator still within the locus?
     * @return True if the iterator has more elements.  False otherwise. 
     */
    public boolean hasNext() {
        return !currentLocus.isPast(completeLocus);    
    }

    /**
     * Get the next single-base locus context bounded by the iterator.
     * @return GenomeLoc representing the next single-base locus context.
     */
    public GenomeLoc next() {
        if( !hasNext() )
            throw new NoSuchElementException("No elements remaining in bounded reference region.");
        GenomeLoc toReturn = currentLocus;
        currentLocus = parser.incPos(currentLocus);
        return toReturn;
    }

    public void remove() {
        throw new UnsupportedOperationException( "ReferenceLocusIterator is read-only" );
    }
}
