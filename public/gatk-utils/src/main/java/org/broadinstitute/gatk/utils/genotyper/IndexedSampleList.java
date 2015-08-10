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

package org.broadinstitute.gatk.utils.genotyper;

import org.broadinstitute.gatk.utils.collections.IndexedSet;

import java.util.Collection;

/**
 * Simple implementation of a sample-list using and indexed-set.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class IndexedSampleList implements SampleList {

    private final IndexedSet<String> samples;

    /**
     * Constructs an empty sample-list.
     */
    public IndexedSampleList() {
        samples = new IndexedSet<>(0);
    }

    /**
     * Constructs a sample-list from a collection of samples.
     *
     * <p>
     *     Repeats in the input collection are ignored (just the first occurrence is kept).
     *     Sample names will be sorted based on the traversal order
     *     of the original collection.
     * </p>
     *
     * @param samples input sample collection.
     *
     * @throws IllegalArgumentException if {@code samples} is {@code null} or it contains {@code nulls}.
     */
    public IndexedSampleList(final Collection<String> samples) {
        this.samples = new IndexedSet<>(samples);
    }

    /**
     * Constructs a sample-list from an array of samples.
     *
     * <p>
     *     Repeats in the input array are ignored (just the first occurrence is kept).
     *     Sample names will be sorted based on the traversal order
     *     of the original array.
     * </p>
     *
     * @param samples input sample array.
     *
     * @throws IllegalArgumentException if {@code samples} is {@code null} or it contains {@code nulls}.
     */
    public IndexedSampleList(final String ... samples) {
        this.samples = new IndexedSet<>(samples);
    }

    @Override
    public int sampleCount() {
        return samples.size();
    }

    @Override
    public int sampleIndex(final String sample) {
        return samples.indexOf(sample);
    }

    @Override
    public String sampleAt(int sampleIndex) {
        return samples.get(sampleIndex);
    }
}
