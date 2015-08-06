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

package org.broadinstitute.gatk.utils.downsampling;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.ArrayList;
import java.util.List;

/**
 * Fractional Downsampler: selects a specified fraction of the reads for inclusion.
 *
 * Since the selection is done randomly, the actual fraction of reads retained may be slightly
 * more or less than the requested fraction, depending on the total number of reads submitted.
 *
 * @author David Roazen
 */
public class FractionalDownsampler<T extends SAMRecord> extends ReadsDownsampler<T> {

    private ArrayList<T> selectedReads;

    private final int cutoffForInclusion;

    private static final int RANDOM_POOL_SIZE = 10000;

    /**
     * Construct a FractionalDownsampler
     *
     * @param fraction Fraction of reads to preserve, between 0.0 (inclusive) and 1.0 (inclusive).
     *                 Actual number of reads preserved may differ randomly.
     */
    public FractionalDownsampler( final double fraction ) {
        if ( fraction < 0.0 || fraction > 1.0 ) {
            throw new ReviewedGATKException("Fraction of reads to include must be between 0.0 and 1.0, inclusive");
        }

        cutoffForInclusion = (int)(fraction * RANDOM_POOL_SIZE);
        clearItems();
        resetStats();
    }

    @Override
    public void submit( final T newRead ) {
        if ( Utils.getRandomGenerator().nextInt(10000) < cutoffForInclusion || doNotDiscardItem(newRead) ) {
            selectedReads.add(newRead);
        }
        else {
            numDiscardedItems++;
        }
    }

    @Override
    public boolean hasFinalizedItems() {
        return selectedReads.size() > 0;
    }

    @Override
    public List<T> consumeFinalizedItems() {
        // pass by reference rather than make a copy, for speed
        List<T> downsampledItems = selectedReads;
        clearItems();
        return downsampledItems;
    }

    @Override
    public boolean hasPendingItems() {
        return false;
    }

    @Override
    public T peekFinalized() {
        return selectedReads.isEmpty() ? null : selectedReads.get(0);
    }

    @Override
    public T peekPending() {
        return null;
    }

    @Override
    public int size() {
        return selectedReads.size();
    }

    @Override
    public void signalEndOfInput() {
        // NO-OP
    }

    @Override
    public void clearItems() {
        selectedReads = new ArrayList<T>();
    }

    @Override
    public boolean requiresCoordinateSortOrder() {
        return false;
    }

    @Override
    public void signalNoMoreReadsBefore( final T read ) {
        // NO-OP
    }
}
