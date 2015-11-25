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

import htsjdk.variant.variantcontext.Allele;

import java.util.AbstractList;
import java.util.List;

/**
 * Utils operations on {@link AlleleList} instances.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class AlleleListUtils {

    @SuppressWarnings("unchecked")
    private static final AlleleList EMPTY_LIST = new AlleleList() {
        @Override
        public int alleleCount() {
            return 0;
        }

        @Override
        public int alleleIndex(final Allele allele) {
            return -1;
        }

        @Override
        public Allele alleleAt(final int index) {
            throw new IllegalArgumentException("allele index is out of range");
        }
    };

    /**
     * Checks whether two allele lists are in fact the same.
     * @param first one list to compare.
     * @param second another list to compare.
     *
     * @throws IllegalArgumentException if if either list is {@code null}.
     *
     * @return {@code true} iff both list are equal.
     */
    public static <A extends Allele> boolean equals(final AlleleList<A> first, final AlleleList<A> second) {
        if (first == null || second == null)
            throw new IllegalArgumentException("no null list allowed");
        final int alleleCount = first.alleleCount();
        if (alleleCount != second.alleleCount())
            return false;

        for (int i = 0; i < alleleCount; i++) {
            final A firstSample = first.alleleAt(i);
            if (firstSample == null)
                throw new IllegalStateException("no null samples allowed in sample-lists: first list at " + i);
            final A secondSample = second.alleleAt(i);
            if (secondSample == null)
                throw new IllegalArgumentException("no null samples allowed in sample-list: second list at " + i);
            if (!firstSample.equals(secondSample))
                return false;
        }

        return true;
    }

    /**
     * Resolves the index of the reference allele in an allele-list.
     *
     * <p>
     *     If there is no reference allele, it returns -1. If there is more than one reference allele,
     *     it returns the first occurrence (lowest index).
     * </p>
     *
     * @param list the search allele-list.
     * @param <A> allele component type.
     *
     * @throws IllegalArgumentException if {@code list} is {@code null}.
     *
     * @return -1 if there is no reference allele, or a values in [0,{@code list.alleleCount()}).
     */
    public static <A extends Allele> int indexOfReference(final AlleleList<A> list) {
        if (list == null)
            throw new IllegalArgumentException("the input list cannot be null");
        final int alleleCount = list.alleleCount();
        for (int i = 0; i < alleleCount; i++)
            if (list.alleleAt(i).isReference())
                return i;
        return -1;
    }


    /**
     * Returns a {@link java.util.List} unmodifiable view of a allele-list
     * @param list the sample-list to wrap.
     *
     * @throws IllegalArgumentException if {@code list} is {@code null}.
     *
     * @return never {@code null}.
     */
    public static <A extends Allele> List<A> asList(final AlleleList<A> list) {
        if (list == null)
            throw new IllegalArgumentException("the list cannot be null");
        return new AsList(list);
    }

    /**
     * Returns an unmodifiable empty allele-list.
     * @param <A> the allele class.
     * @return never {@code null}.
     */
    @SuppressWarnings("unchecked")
    public static final <A extends Allele> AlleleList<A> emptyList() {
        return EMPTY_LIST;
    }

    /**
     * Simple list view of a sample-list.
     */
    private static class AsList<A extends Allele> extends AbstractList<A> {

        private final AlleleList<A> list;

        private AsList(final AlleleList<A> list) {
            this.list = list;

        }

        @Override
        public A get(int index) {
            return list.alleleAt(index);
        }

        @Override
        public int size() {
            return list.alleleCount();
        }
    }


    /**
     * Returns a permutation between two allele lists.
     * @param original the original allele list.
     * @param target the target allele list.
     * @param <A> the allele type.
     *
     * @throws IllegalArgumentException if {@code original} or {@code target} is {@code null}, or
     * elements in {@code target} is not contained in {@code original}
     *
     * @return never {@code null}
     */
    public static <A extends Allele> AlleleListPermutation<A> permutation(final AlleleList<A> original, final AlleleList<A> target) {
        if (equals(original,target))
            return new NonPermutation<>(original);
        else
            return new ActualPermutation<>(original,target);
    }

    private static class NonPermutation<A extends Allele> implements AlleleListPermutation<A> {

        private final AlleleList<A> list;

        public NonPermutation(final AlleleList<A> original) {
            list = original;
        }

        @Override
        public boolean isPartial() {
            return false;
        }

        @Override
        public boolean isNonPermuted() {
            return true;
        }

        @Override
        public int toIndex(int fromIndex) {
            return fromIndex;
        }

        @Override
        public int fromIndex(int toIndex) {
            return toIndex;
        }

        @Override
        public int fromSize() {
            return list.alleleCount();
        }

        @Override
        public int toSize() {
            return list.alleleCount();
        }

        @Override
        public List<A> fromList() {
            return asList(list);
        }

        @Override
        public java.util.List<A> toList() {
            return asList(list);
        }


        @Override
        public int alleleCount() {
            return list.alleleCount();
        }

        @Override
        public int alleleIndex(final A allele) {
            return list.alleleIndex(allele);
        }

        @Override
        public A alleleAt(final int index) {
            return list.alleleAt(index);
        }
    }

    private static class ActualPermutation<A extends Allele> implements AlleleListPermutation<A> {

        private final AlleleList<A> from;

        private final AlleleList<A> to;

        private final int[] fromIndex;

        private final boolean nonPermuted;

        private final boolean isPartial;

        private ActualPermutation(final AlleleList<A> original, final AlleleList<A> target) {
            this.from = original;
            this.to = target;
            final int toSize = target.alleleCount();
            final int fromSize = original.alleleCount();
            if (fromSize < toSize)
                throw new IllegalArgumentException("target allele list is not a permutation of the original allele list");

            fromIndex = new int[toSize];
            boolean nonPermuted = fromSize == toSize;
            this.isPartial = !nonPermuted;
            for (int i = 0; i < toSize; i++) {
                final int originalIndex = original.alleleIndex(target.alleleAt(i));
                if (originalIndex < 0)
                    throw new IllegalArgumentException("target allele list is not a permutation of the original allele list");
                fromIndex[i] = originalIndex;
                nonPermuted &= originalIndex == i;
            }

            this.nonPermuted = nonPermuted;
        }

        @Override
        public boolean isPartial() {
            return isPartial;
        }

        @Override
        public boolean isNonPermuted() {
            return nonPermuted;
        }

        @Override
        public int toIndex(int fromIndex) {
            return to.alleleIndex(from.alleleAt(fromIndex));
        }

        @Override
        public int fromIndex(int toIndex) {
            return fromIndex[toIndex];
        }

        @Override
        public int fromSize() {
            return from.alleleCount();
        }

        @Override
        public int toSize() {
            return to.alleleCount();
        }

        @Override
        public List<A> fromList() {
            return asList(from);
        }

        @Override
        public List<A> toList() {
            return asList(to);
        }

        @Override
        public int alleleCount() {
            return to.alleleCount();
        }

        @Override
        public int alleleIndex(final A allele) {
            return to.alleleIndex(allele);
        }

        @Override
        public A alleleAt(final int index) {
            return to.alleleAt(index);
        }
    }
}
