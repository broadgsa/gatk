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

import java.util.*;

/**
 * Some utility operations on sample lists.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class SampleListUtils {

    private static final SampleList EMPTY_LIST = new SampleList() {

        @Override
        public int sampleCount() {
            return 0;
        }

        @Override
        public int sampleIndex(String sample) {
            return -1;
        }

        @Override
        public String sampleAt(final int sampleIndex) {
            throw new IllegalArgumentException("index is out of valid range");
        }
    };

    /**
     * Empty list.
     *
     * @return never {@code null}
     */
    public static SampleList emptyList() {
        return EMPTY_LIST;
    }

    /**
     * Checks whether two sample lists are in fact the same.
     * @param first one list to compare.
     * @param second another list to compare.
     *
     * @throws IllegalArgumentException if if either list is {@code null}.
     *
     * @return {@code true} iff both list are equal.
     */
    public static boolean equals(final SampleList first, final SampleList second) {
        if (first == null || second == null)
            throw new IllegalArgumentException("no null list allowed");
        final int sampleCount = first.sampleCount();
        if (sampleCount != second.sampleCount())
            return false;

        for (int i = 0; i < sampleCount; i++) {
            final String firstSample = first.sampleAt(i);
            if (firstSample == null)
                throw new IllegalStateException("no null samples allowed in sample-lists: first list at " + i);
            final String secondSample = second.sampleAt(i);
            if (secondSample == null)
                throw new IllegalArgumentException("no null samples allowed in sample-list: second list at " + i);
            if (!firstSample.equals(secondSample))
                return false;
        }
        return true;
    }

    /**
     * Returns a {@link List} unmodifiable view of a sample-list
     * @param list the sample-list to wrap.
     *
     * @throws IllegalArgumentException if {@code list} is {@code null}.
     *
     * @return never {@code null}.
     */
    public static List<String> asList(final SampleList list) {
        if (list == null)
            throw new IllegalArgumentException("the list cannot be null");
        return new AsList(list);
    }

    /**
     * Returns a {@link Set} unmodifiable view of the sample-list
     *
     * @param list the sample-list to wrap.
     *
     * @throws IllegalArgumentException if {@code list} is {@code null}
     */
    public static Set<String> asSet(final SampleList list) {
        if (list == null)
            throw new IllegalArgumentException("the list cannot be null");
        return new AsSet(list);
    }

    /**
     * Creates a list with a single sample.
     *
     * @param sampleName the sample name.
     * @return never {@code sampleName}
     */
    public static SampleList singletonList(final String sampleName) {
        if (sampleName == null)
            throw new IllegalArgumentException("the sample name cannot be null");
        return new SampleList() {

            @Override
            public int sampleCount() {
                return 1;
            }

            @Override
            public int sampleIndex(final String sample) {
                return sampleName.equals(sample) ? 0 : -1;
            }

            @Override
            public String sampleAt(int sampleIndex) {
                if (sampleIndex == 0)
                    return sampleName;
                throw new IllegalArgumentException("index is out of bounds");
            }
        };
    }

    /**
     * Simple list view of a sample-list.
     */
    private static class AsList extends AbstractList<String> {

        private final SampleList list;

        private AsList(final SampleList list) {
            this.list = list;

        }

        @Override
        public String get(int index) {
            return list.sampleAt(index);
        }

        @Override
        public int size() {
            return list.sampleCount();
        }
    }

    /**
     * Simple set view of a sample-list
     */
    private static class AsSet extends AbstractSet<String> {

        private final SampleList list;

        private AsSet(final SampleList list) {
            this.list = list;

        }

        @Override
        public Iterator<String> iterator() {
            return new Iterator<String>() {
                private int index = 0;

                @Override
                public boolean hasNext() {
                    return index < list.sampleCount();
                }

                @Override
                public String next() {
                    if (index >= list.sampleCount())
                        throw new NoSuchElementException("iterating beyond sample list end");
                    return list.sampleAt(index++);
                }

                @Override
                public void remove() {
                    throw new UnsupportedOperationException("unsupported operation exception");
                }
            };
        }

        @Override
        public int size() {
            return list.sampleCount();
        }

        @Override
        public boolean contains(final Object obj) {
            if (obj == null)
                return false;
            else if (obj instanceof String)
                return list.sampleIndex(((String)obj)) >= 0;
            else
                return false;
        }
    }
}
