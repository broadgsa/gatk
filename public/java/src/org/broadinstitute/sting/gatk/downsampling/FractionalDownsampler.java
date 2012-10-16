/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.downsampling;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Fractional Downsampler: selects a specified fraction of the reads for inclusion.
 *
 * Since the selection is done randomly, the actual fraction of reads retained may be slightly
 * more or less than the requested fraction, depending on the total number of reads submitted.
 *
 * @author David Roazen
 */
public class FractionalDownsampler<T extends SAMRecord> implements ReadsDownsampler<T> {

    private ArrayList<T> selectedReads;

    private int cutoffForInclusion;

    private int numDiscardedItems;

    private static final int RANDOM_POOL_SIZE = 10000;

    /**
     * Construct a FractionalDownsampler
     *
     * @param fraction Fraction of reads to preserve, between 0.0 (inclusive) and 1.0 (inclusive).
     *                 Actual number of reads preserved may differ randomly.
     */
    public FractionalDownsampler( double fraction ) {
        if ( fraction < 0.0 || fraction > 1.0 ) {
            throw new ReviewedStingException("Fraction of reads to include must be between 0.0 and 1.0, inclusive");
        }

        cutoffForInclusion = (int)(fraction * RANDOM_POOL_SIZE);
        clear();
        reset();
    }

    public void submit( T newRead ) {
        if ( GenomeAnalysisEngine.getRandomGenerator().nextInt(10000) < cutoffForInclusion ) {
            selectedReads.add(newRead);
        }
        else {
            numDiscardedItems++;
        }
    }

    public void submit( Collection<T> newReads ) {
        for ( T read : newReads ) {
            submit(read);
        }
    }

    public boolean hasFinalizedItems() {
        return selectedReads.size() > 0;
    }

    public List<T> consumeFinalizedItems() {
        // pass by reference rather than make a copy, for speed
        List<T> downsampledItems = selectedReads;
        clear();
        return downsampledItems;
    }

    public boolean hasPendingItems() {
        return false;
    }

    public T peekFinalized() {
        return selectedReads.isEmpty() ? null : selectedReads.get(0);
    }

    public T peekPending() {
        return null;
    }

    public int getNumberOfDiscardedItems() {
        return numDiscardedItems;
    }

    public void signalEndOfInput() {
        // NO-OP
    }

    public void clear() {
        selectedReads = new ArrayList<T>();
    }

    public void reset() {
        numDiscardedItems = 0;
    }

    public boolean requiresCoordinateSortOrder() {
        return false;
    }

    public void signalNoMoreReadsBefore( T read ) {
        // NO-OP
    }
}
