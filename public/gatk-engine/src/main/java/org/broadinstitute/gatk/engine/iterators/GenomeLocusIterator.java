/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.iterators;

import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;

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
