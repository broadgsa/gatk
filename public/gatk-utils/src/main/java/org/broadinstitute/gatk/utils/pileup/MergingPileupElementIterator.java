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

package org.broadinstitute.gatk.utils.pileup;

import htsjdk.samtools.util.PeekableIterator;

import java.util.Comparator;
import java.util.Iterator;
import java.util.PriorityQueue;

/**
 * Merges multiple pileups broken down by sample.
 *
 * @author mhanna
 * @version 0.1
 */
class MergingPileupElementIterator<PE extends PileupElement> implements Iterator<PE> {
    private final PriorityQueue<PeekableIterator<PE>> perSampleIterators;

    public MergingPileupElementIterator(PerSamplePileupElementTracker<PE> tracker) {
        perSampleIterators = new PriorityQueue<PeekableIterator<PE>>(Math.max(1,tracker.getSamples().size()),new PileupElementIteratorComparator());
        for(final String sample: tracker.getSamples()) {
            PileupElementTracker<PE> trackerPerSample = tracker.getElements(sample);
            if(trackerPerSample.size() != 0)
                perSampleIterators.add(new PeekableIterator<PE>(trackerPerSample.iterator()));
        }
    }

    public boolean hasNext() {
        return !perSampleIterators.isEmpty();
    }

    public PE next() {
        PeekableIterator<PE> currentIterator = perSampleIterators.remove();
        PE current = currentIterator.next();
        if(currentIterator.hasNext())
            perSampleIterators.add(currentIterator);
        return current;
    }

    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a merging iterator.");
    }

    /**
     * Compares two peekable iterators consisting of pileup elements.
     */
    private class PileupElementIteratorComparator implements Comparator<PeekableIterator<PE>> {
        public int compare(PeekableIterator<PE> lhs, PeekableIterator<PE> rhs) {
            return rhs.peek().getOffset() - lhs.peek().getOffset();
        }
    }
}
