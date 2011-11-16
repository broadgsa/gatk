/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.utils.pileup;

import java.util.*;

/**
 * Javadoc goes here.
 *
 * @author mhanna
 * @version 0.1
 */
abstract class PileupElementTracker<PE extends PileupElement> implements Iterable<PE> {
    public abstract int size();
}

class UnifiedPileupElementTracker<PE extends PileupElement> extends PileupElementTracker<PE> {
    private final List<PE> pileup;

    public UnifiedPileupElementTracker() { pileup = new LinkedList<PE>(); }
    public UnifiedPileupElementTracker(List<PE> pileup) { this.pileup = pileup; }

    public void add(PE element) {
        pileup.add(element);
    }

    public PE get(int index) {
        return pileup.get(index);
    }

    public int size() {
        return pileup.size();
    }

    public Iterator<PE> iterator() { return pileup.iterator(); }
}

class PerSamplePileupElementTracker<PE extends PileupElement> extends PileupElementTracker<PE> {
    private final Map<String,PileupElementTracker<PE>> pileup;
    private int size = 0;

    public PerSamplePileupElementTracker() {
        pileup = new HashMap<String,PileupElementTracker<PE>>();
    }

    /**
     * Gets a list of all the samples stored in this pileup.
     * @return List of samples in this pileup.
     */
    public Collection<String> getSamples() {
        return pileup.keySet();
    }

    public PileupElementTracker<PE> getElements(final String sample) {
        return pileup.get(sample);
    }

    public PileupElementTracker<PE> getElements(final Collection<String> selectSampleNames) {
        PerSamplePileupElementTracker<PE> result = new PerSamplePileupElementTracker<PE>();
        for (final String sample :  selectSampleNames) {
            result.addElements(sample, pileup.get(sample));
        }
        return result;
    }

    public void addElements(final String sample, PileupElementTracker<PE> elements) {
        pileup.put(sample,elements);
        size += elements.size();
    }

    public Iterator<PE> iterator() { return new MergingPileupElementIterator<PE>(this); }

    public int size() {
        return size;
    }
}