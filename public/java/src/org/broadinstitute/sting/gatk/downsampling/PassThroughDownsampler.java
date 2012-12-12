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

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Pass-Through Downsampler: Implementation of the ReadsDownsampler interface that does no
 * downsampling whatsoever, and instead simply "passes-through" all the reads it's given.
 * Useful for situations where you want to disable downsampling, but still need to use
 * the downsampler interface.
 *
 * @author David Roazen
 */
public class PassThroughDownsampler<T extends SAMRecord> implements ReadsDownsampler<T> {

    private ArrayList<T> selectedReads;

    public PassThroughDownsampler() {
        clear();
    }

    public void submit( T newRead ) {
        // All reads pass-through, no reads get downsampled
        selectedReads.add(newRead);
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
        return 0;
    }

    public void signalEndOfInput() {
        // NO-OP
    }

    public void clear() {
        selectedReads = new ArrayList<T>();
    }

    public void reset() {
        // NO-OP
    }

    public boolean requiresCoordinateSortOrder() {
        return false;
    }

    public void signalNoMoreReadsBefore( T read ) {
        // NO-OP
    }
}
