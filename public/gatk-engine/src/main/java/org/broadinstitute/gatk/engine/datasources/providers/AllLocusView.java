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

import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.engine.iterators.GenomeLocusIterator;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.Collections;
import java.util.List;
import java.util.NoSuchElementException;
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
     * Signal not to advance the iterator because we're currently sitting at the next element.
     */
    private boolean atNextElement = false;

    /**
     * Create a new queue of locus contexts.
     *
     * @param provider
     */
    public AllLocusView(LocusShardDataProvider provider) {
        super(provider);
        // Seed the state tracking members with the first possible seek position and the first possible locus context.
        locusIterator = new GenomeLocusIterator(genomeLocParser, provider.getLocus());
    }

    public boolean hasNext() {
        advance();
        return nextPosition != null;
    }

    public AlignmentContext next() {
        advance();

        if (nextPosition == null)
            throw new NoSuchElementException("No next is available in the all locus view");

        // Flag to the iterator that no data is waiting in the queue to be processed.
        atNextElement = false;

        AlignmentContext currentLocus;

        // If actual data is present, return it.  Otherwise, return empty data.
        if (nextLocus != null && nextLocus.getLocation().equals(nextPosition))
            currentLocus = nextLocus;
        else
            currentLocus = createEmptyLocus(nextPosition);

        return currentLocus;
    }

    private void advance() {
        // Already at the next element?  Don't move forward.
        if (atNextElement)
            return;

        // Out of elements?
        if (nextPosition == null && !locusIterator.hasNext())
            return;

        // If nextLocus has been consumed, clear it out to make room for the next incoming locus.
        if (nextPosition != null && nextLocus != null && !nextLocus.getLocation().isPast(nextPosition)) {
            nextLocus = null;

            // Determine the next locus. The trick is that we may have more than one alignment context at the same
            // reference position (regular base pileup, then extended pileup). If next alignment context (that we just pre-read)
            // is still at the current position, we do not increment current position and wait for next call to next() to return
            // that context. If we know that next context is past the current position, we are done with current
            // position
            if (hasNextLocus()) {
                nextLocus = nextLocus();
                if (nextPosition.equals(nextLocus.getLocation())) {
                    atNextElement = true;
                    return;
                }
            }
        }

        // No elements left in queue?  Clear out the position state tracker and return.
        if (!locusIterator.hasNext()) {
            nextPosition = null;
            return;
        }

        // Actually fill the next position.
        nextPosition = locusIterator.next();
        atNextElement = true;

        // Crank the iterator to (if possible) or past the next context.  Be careful not to hold a reference to nextLocus
        // while using the hasNextLocus() / nextLocus() machinery; this will cause us to use more memory than is optimal. 
        while (nextLocus == null || nextLocus.getLocation().isBefore(nextPosition)) {
            nextLocus = null;
            if (!hasNextLocus())
                break;
            nextLocus = nextLocus();
        }
    }

    /**
     * Creates a blank locus context at the specified location.
     *
     * @param site Site at which to create the blank locus context.
     * @return empty context.
     */
    private final static List<GATKSAMRecord> EMPTY_PILEUP_READS = Collections.emptyList();
    private final static List<Integer> EMPTY_PILEUP_OFFSETS = Collections.emptyList();
    private final static List<Boolean> EMPTY_DELETION_STATUS = Collections.emptyList();

    private AlignmentContext createEmptyLocus(GenomeLoc site) {
        return new AlignmentContext(site, new ReadBackedPileupImpl(site, EMPTY_PILEUP_READS, EMPTY_PILEUP_OFFSETS));
    }
}
