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

package org.broadinstitute.gatk.engine.datasources.providers;

import org.broadinstitute.gatk.engine.ReadProperties;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.locusiterator.LocusIterator;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.locusiterator.LocusIteratorByState;

import java.util.Arrays;
import java.util.Collection;
import java.util.NoSuchElementException;

/**
 * User: hanna
 * Date: May 13, 2009
 * Time: 3:30:16 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * The two goals of the LocusView are as follows:
 * 1) To provide a 'trigger track' iteration interface so that TraverseLoci can easily switch
 *    between iterating over all bases in a region, only covered bases in a region covered by
 *    reads, only bases in a region covered by RODs, or any other sort of trigger track
 *    implementation one can think of.
 * 2) To manage the copious number of iterators that have to be jointly pulled through the
 *    genome to make a locus traversal function.
 */
public abstract class LocusView extends LocusIterator implements View {
    /**
     * The locus bounding this view.
     */
    protected GenomeLoc locus;

    /**
     * The GenomeLocParser, used to create new genome locs.
     */
    protected GenomeLocParser genomeLocParser;

    /**
     * Source info for this view.  Informs the class about downsampling requirements.
     */
    private ReadProperties sourceInfo;

    /**
     * The actual locus context iterator.
     */
    private LocusIterator loci;

    /**
     * The next locus context from the iterator.  Lazy loaded: if nextLocus is null and advance() doesn't
     * populate it, the iterator is exhausted.  If populated, this is the value that should be returned by
     * next(). 
     */
    private AlignmentContext nextLocus = null;

    public LocusView(LocusShardDataProvider provider) {
        this.locus = provider.getLocus();
        
        this.sourceInfo = provider.getSourceInfo();
        this.genomeLocParser = provider.getGenomeLocParser();
        this.loci = provider.getLocusIterator();

        advance();

        provider.register(this);
    }

    /**
     * Only one view of the locus is supported at any given time.
     * @return A list consisting of all other locus views.
     */
    public Collection<Class<? extends View>> getConflictingViews() {
        return Arrays.<Class<? extends View>>asList(LocusView.class,ReadView.class);
    }

    /**
     * Close this view.
     */
    public void close() {
        // Set everything to null with the hope of failing fast.
        locus = null;
        sourceInfo = null;
        loci = null;

        super.close();
    }

    /**
     * Is there another covered locus context bounded by this view.
     * @return True if another covered locus context exists.  False otherwise.
     */
    public abstract boolean hasNext();

    /**
     * Returns the next covered locus context in the shard.
     * @return Next covered locus context in the shard.
     * @throw NoSuchElementException if no such element exists.
     */
    public abstract AlignmentContext next();

    /**
     * Unsupported.
     * @throw UnsupportedOperationException always.
     */
    public void remove() {
        throw new UnsupportedOperationException("Unable to remove elements from this queue.");
    }

    /**
     * Is there another locus context bounded by this shard.
     * @return True if another locus context is bounded by this shard.
     */
    protected boolean hasNextLocus() {
        advance();
        return nextLocus != null;
    }

    /**
     * Get the next locus context bounded by this shard.
     * @return Next locus context bounded by this shard.
     * @throw NoSuchElementException if the next element is missing.
     */
    protected AlignmentContext nextLocus() {
        advance();
        if(nextLocus == null)
            throw new NoSuchElementException("No more elements remain in locus context queue.");

        // Cache the current and apply filtering.
        AlignmentContext current = nextLocus;

        // Indicate that the next operation will need to advance.
        nextLocus = null;
        
        return current;
    }

    /**
     * Seed the nextLocus variable with the contents of the next locus (if one exists).
     */
    private void advance() {
        // Already an unclaimed locus present
        if(nextLocus != null)
            return;

        //System.out.printf("loci is %s%n", loci);
        if( !loci.hasNext() ) {
            nextLocus = null;
            return;
        }

        nextLocus = loci.next();

        // If the location of this shard is available, trim the data stream to match the shard.
        // TODO: Much of this functionality is being replaced by the WindowMaker.
        if(locus != null) {
            // Iterate through any elements not contained within this shard.
            while( nextLocus != null && !isContainedInShard(nextLocus.getLocation()) && loci.hasNext() )
                nextLocus = loci.next();

            // If nothing in the shard was found, indicate that by setting nextLocus to null.
            if( nextLocus != null && !isContainedInShard(nextLocus.getLocation()) )
                nextLocus = null;
        }
    }

    /**
     * Is this location contained in the given shard.
     * @param location Location to check.
     * @return True if the given location is contained within the shard.  False otherwise.
     */
    private boolean isContainedInShard(GenomeLoc location) {
        return locus.containsP(location);
    }

    /**
     * {@inheritDoc}
     *
     * Since this class has an actual LIBS, so this function will never throw an exception
     *
     * @return the LocusIteratorByState used by this view to get pileups
     */
    @Override
    public LocusIteratorByState getLIBS() {
        return loci.getLIBS();
    }
}
